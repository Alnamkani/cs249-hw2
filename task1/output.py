def output():
    metrics_dict = {
    "Sequence length": "Total assembly length",
    "Number of contigs": "Total number of contigs",
    "GC content (%)": "GC percentage",
    "Genome fraction (%)": "% of reference covered by assembly",
    "Duplication ratio": "Average number of times a reference base is covered",
    "Largest contig": "Length of the largest contig",
    "N50": "Length where 50% of assembly is in contigs of this size or larger",
    "N90": "Length where 90% of assembly is in contigs of this size or larger",
    "L50": "Number of contigs to reach 50% of total assembly length",
    "Misassemblies": "Number of positions with breakpoints relative to reference",
    "Mismatches per 100 kbp": "Number of mismatches per 100,000 bases",
    "Indels per 100 kbp": "Number of insertions/deletions per 100,000 bases"
}