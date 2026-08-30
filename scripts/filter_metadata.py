import pandas as pd

counts_header = pd.read_csv(
    snakemake.input["counts"], sep="\t", nrows=0
).columns.tolist()
included_samples = set(counts_header) - {"gene_id"}

metadata = pd.read_excel(snakemake.input["metadata"])
filtered = metadata[metadata["Sample_ID"].astype(str).isin(included_samples)]

if filtered.empty:
    raise ValueError(
        "filter_metadata produced zero rows -- check that Sample_ID values in "
        "the metadata file match the column headers in gene_raw_counts.txt."
    )

filtered.to_excel(snakemake.output[0], index=False)
