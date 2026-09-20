#!/usr/bin/env python3
"""
Critic / verifier layer
"""

import re
from typing import Optional
from dataclasses import dataclass


@dataclass
class VerificationResult:
    gene_id: str
    supported: bool
    reason: str
    check: str
    evidence: Optional[str] = None

    def as_dict(self) -> dict:
        return {
            "gene_id": self.gene_id,
            "supported": self.supported,
            "reason": self.reason,
            "check": self.check,
            "evidence": self.evidence,
        }


def _check_claim_not_empty(
    gene_id: str, claimed_function: str
) -> Optional[VerificationResult]:
    if not claimed_function or not claimed_function.strip():
        return VerificationResult(
            gene_id=gene_id,
            supported=False,
            reason="Claimed function is empty -- nothing to verify.",
            check="claim_nonempty",
        )
    return None


def _check_gene_id_was_retrieved(
    gene_id: str, annotation_record: dict
) -> Optional[VerificationResult]:
    if not gene_id or not gene_id.strip():
        return VerificationResult(
            gene_id=gene_id,
            supported=False,
            reason="No gene_id given -- claim isn't attributable to anything.",
            check="gene_id_retrieved",
        )
    if not annotation_record.get("found"):
        return VerificationResult(
            gene_id=gene_id,
            supported=False,
            reason=f"get_gene_annotation found nothing for {gene_id!r} -- claim has nothing to be grounded in.",
            check="gene_id_retrieved",
        )
    if annotation_record.get("gene_id") != gene_id:
        return VerificationResult(
            gene_id=gene_id,
            supported=False,
            reason=f"Retrieved record is for {annotation_record.get('gene_id')!r}, not {gene_id!r} -- claim may be about the wrong gene.",
            check="gene_id_retrieved",
        )
    return None


_MECHANICAL_CHECKS = (
    _check_claim_not_empty,
    _check_gene_id_was_retrieved,
)
_BI_ENCODER_MODEL_NAME = "sentence-transformers/all-MiniLM-L6-v2"
_NLI_MODEL_NAME = "cross-encoder/nli-deberta-v3-small"

# Starting points only -- calibrate against evals/cases/ once you can run it.
_COSINE_PREFILTER_THRESHOLD = 0.35
_ENTAILMENT_THRESHOLD = 0.5

_HEDGE_MARKERS = (
    "putative",
    "predicted",
    "inferred",
    "candidate",
    "possible",
    "probable",
    "likely",
    "hypothetical",
    "uncharacterized",
    "unconfirmed",
    "proposed",
    "may be",
    "believed to",
    "thought to",
    "suggested to",
)

_bi_encoder = None
_nli_model = None
_nli_id2label = None


def _get_bi_encoder():
    global _bi_encoder
    if _bi_encoder is None:
        from sentence_transformers import SentenceTransformer

        _bi_encoder = SentenceTransformer(_BI_ENCODER_MODEL_NAME)
    return _bi_encoder


def _get_nli_model():
    global _nli_model, _nli_id2label
    if _nli_model is None:
        from sentence_transformers import CrossEncoder

        _nli_model = CrossEncoder(_NLI_MODEL_NAME)
        _nli_id2label = {
            int(k): v.lower() for k, v in _nli_model.model.config.id2label.items()
        }
    return _nli_model, _nli_id2label


def _split_sentences(text: str) -> list:
    text = text.strip()
    protected = re.sub(r"\b([A-Z])\.\s+(?=[a-z])", r"\1<DOT> ", text)
    sentences = re.split(r"(?<=[.!?])\s+", protected)
    return [s.replace("<DOT>", ".").strip() for s in sentences if s.strip()]


def _cosine_similarity(a, b) -> float:
    import numpy as np

    a = np.asarray(a)
    b = np.asarray(b)
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom == 0.0:
        return 0.0
    return float(np.dot(a, b) / denom)


def _most_similar_sentence(claim: str, annotation_sentences: list):
    model = _get_bi_encoder()
    candidates = list(dict.fromkeys(annotation_sentences + [" ".join(annotation_sentences)]))
    claim_vec = model.encode(claim)
    sentence_vecs = model.encode(candidates)
    sims = [_cosine_similarity(claim_vec, v) for v in sentence_vecs]
    best_idx = max(range(len(sims)), key=lambda i: sims[i])
    return candidates[best_idx], sims[best_idx]


def _hedges_present(text: str) -> list:
    text_lower = text.lower()
    return [m for m in _HEDGE_MARKERS if re.search(rf"\b{re.escape(m)}\b", text_lower)]

# Three-pass design, in order:
#   1. Cosine similarity (bi-encoder) as a cheap topical prefilter -- rules
#      out claims that aren't even about the same subject as the retrieved
#      text.
#   2. Cross-encoder NLI checks whether the claim is actually entailed by
#      the closest matching sentence, not just topically related to it.
#   3. An explicit hedge-word check, because NLI alone scores "putative
#      adhesin" -> "adhesin" as entailment (dropping a hedge doesn't
#      technically contradict the source) -- overclaiming needs its own
#      rule, entailment won't catch it.

def _check_semantic_support(
    gene_id: str, claimed_function: str, retrieved_annotation: str
) -> VerificationResult:
    sentences = _split_sentences(retrieved_annotation) or [retrieved_annotation]
    best_sentence, cosine_score = _most_similar_sentence(claimed_function, sentences)

    if cosine_score < _COSINE_PREFILTER_THRESHOLD:
        return VerificationResult(
            gene_id=gene_id,
            supported=False,
            reason=(
                f"Closest match in the annotation has cosine similarity "
                f"{cosine_score:.2f}, below the {_COSINE_PREFILTER_THRESHOLD} "
                f"topical pre-filter -- claim doesn't appear to be about "
                f"anything the annotation actually says. Closest sentence: "
                f"{best_sentence!r}"
            ),
            check="semantic_support_cosine_prefilter",
            evidence=best_sentence,
        )

    nli_model, id2label = _get_nli_model()
    raw_scores = nli_model.predict(
        [(best_sentence, claimed_function)], apply_softmax=True
    )[0]
    label_scores = {id2label[i]: float(raw_scores[i]) for i in range(len(raw_scores))}
    top_label = max(label_scores, key=label_scores.get)
    top_score = label_scores[top_label]

    if top_label != "entailment" or top_score < _ENTAILMENT_THRESHOLD:
        return VerificationResult(
            gene_id=gene_id,
            supported=False,
            reason=(
                f"NLI model's top label was {top_label!r} ({top_score:.2f} "
                f"confidence) against the closest annotation sentence, not "
                f"a confident entailment. Compared against: {best_sentence!r}"
            ),
            check="semantic_support_nli",
            evidence=best_sentence,
        )

    hedges_in_annotation = _hedges_present(best_sentence)
    hedges_in_claim = _hedges_present(claimed_function)
    if hedges_in_annotation and not hedges_in_claim:
        return VerificationResult(
            gene_id=gene_id,
            supported=False,
            reason=(
                f"Annotation hedges this ({', '.join(hedges_in_annotation)} "
                f"present in {best_sentence!r}) but the claim states it with "
                f"no corresponding hedge -- overclaiming certainty the source "
                f"doesn't support. NLI scored this as entailment ({top_score:.2f}) "
                f"because dropping a hedge doesn't strictly contradict the "
                f"source, which is exactly why this needs an explicit check."
            ),
            check="semantic_support_hedge_dropped",
            evidence=best_sentence,
        )

    return VerificationResult(
        gene_id=gene_id,
        supported=True,
        reason=(
            f"NLI model scored this claim as entailment ({top_score:.2f} "
            f"confidence) against the closest annotation sentence, with no "
            f"hedge dropped. Matched: {best_sentence!r}"
        ),
        check="semantic_support_nli",
        evidence=best_sentence,
    )


def verify_claim_against_annotation(
    gene_id: str, claimed_function: str, annotation_record: dict
) -> VerificationResult:
    for check in _MECHANICAL_CHECKS:
        result = check(
            gene_id,
            (
                claimed_function
                if check is _check_claim_not_empty
                else annotation_record
            ),
        )
        if result is not None:
            return result

    return _check_semantic_support(
        gene_id, claimed_function, annotation_record.get("description", "")
    )


def verify_summary(claims: list) -> list:
    return [
        verify_claim_against_annotation(
            c["gene_id"], c["claimed_function"], c["annotation_record"]
        )
        for c in claims
    ]
