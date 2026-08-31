import pandas as pd

counts_header = pd.read_csv(
    snakemake.input["counts"], sep="\t", nrows=0
).columns.tolist()
included_run_accessions = set(counts_header) - {"gene_id"}

metadata = pd.read_excel(snakemake.input["metadata"])
metadata["run_accession"] = metadata["FASTQ_file"].astype(str).str.replace(
    r"\.fastq(\.gz)?$", "", regex=True
)
filtered = metadata[metadata["run_accession"].isin(included_run_accessions)]

if filtered.empty:
    raise ValueError(
        "filter_metadata produced zero rows -- check that FASTQ_file-derived "
        "run accessions in the metadata file match the column headers in "
        "gene_raw_counts.txt."
    )

filtered.to_excel(snakemake.output[0], index=False)
