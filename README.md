# Pseudogenes-Metabolic-Integration

Pipeline to detect candidate pseudogenes in MAGs (metagenome-assembled genomes) and integrate pseudogene predictions with metabolic/functional annotations.

## Repository structure

```
input/
  My_reference_cds.fna
  My_MAG_genomic.fa
  prokka_out/
    My_MAG.ffn
    My_MAG.faa
    emapper/
      My_MAG.emapper.annotations

tools/
  detect_stops_compare_length.py
  remap_coverage_variants.sh
  compare_pseudogenes_with_reference.sh
  incorporate_annotations.py
  compare_manual_with_pseudofinder.py

Pseudofinder/
  My_MAG_pseudos.gff

metabolism/
  pseudogenes_metabolism_integration.ipynb
```

## Files of interest

- `requirements.txt` — Python dependencies (Biopython, pandas, numpy)
- `pseudogenes_metabolism_workflow.txt` — More detailed workflow notes
- `metabolism/pseudogenes_metabolism_integration.ipynb` — Jupyter Notebook for data analysis

---

## Quick example / workflow

1) Annotate with Prokka to get CDS and proteins:

```bash
prokka --outdir prokka_out --prefix My_MAG --addgenes --addmrna --locustag MYMAG \
  --compliant --centre mRs --genus Genus --metagenome --evalue 1e-05 My_MAG_genomic.fa
```

2) Annotate proteins with EggNOG-mapper (emapper):

```bash
emapper.py --cpu 32 \
  -m diamond \
  -i My_MAG_proteins.faa \
  -o My_MAG \
  --output_dir ./emapper \
  --report_orthologs --decorate_gff yes
```

3) Detect internal stop codons and compare lengths (MAG vs reference proteins)

Usage/help:

```bash
python3 tools/detect_stops_compare_length.py -h
usage: detect_stops_compare_length.py [-h] --cds_fasta CDS_FASTA --mag_faa MAG_FAA --ref_faa REF_FAA --prefix PREFIX

Detect internal stop codons, compare lengths vs reference, and summarize likely pseudogenes.

options:
  -h, --help            show this help message and exit
  --cds_fasta CDS_FASTA
                        Path to the nucleotide CDS FASTA file (e.g., Prokka .ffn).
  --mag_faa MAG_FAA     Path to the MAG proteins FASTA file.
  --ref_faa REF_FAA     Path to the reference proteins FASTA file.
  --prefix PREFIX       Prefix for output files.
```

4) Validate pseudogene candidates using read mapping and variant calls

- Remap the MAG CDS (Prokka `.ffn`) against the source metagenome using Bowtie2 and samtools
- Call variants with `bcftools`
- Get read coverage and filter variants by quality

Usage/help:

```bash
./tools/remap_coverage_variants.sh -h
Usage: remap_coverage_variants.sh --ffn_file <path> --metagenome_file <path> --prefix <prefix>

Arguments:
  --ffn_file           Path to the nucleotide CDS FASTA file (e.g., Prokka .ffn).
  --metagenome_file    Path to the metagenome file.
  --prefix             Prefix for output files.
```

5) Analyze alignments between pseudogene candidates and the reference genome using `blastn` (genomic sequences)

Look for:
- Truncated proteins (alignment length < 70–80% of reference)
- Frameshifts or internal stop codons (high mismatch or gapopen values)
- Percent identity (high identity suggests conserved regions)
- Parse BLAST output to identify pseudogenes with disruptions (truncated alignments, frameshifts, etc.)

Usage/help:

```bash
./tools/compare_pseudogenes_with_reference.sh -h
Usage: compare_pseudogenes_with_reference.sh --ffn_file <path> --reference_file <path> --stops_length_file <path> --prefix <prefix>

Arguments:
  --ffn_file           Path to the nucleotide CDS FASTA file (e.g., Prokka .ffn).
  --reference_file     Path to the reference genome file.
  --stops_length_file  Path to the mRsS25-Chr1-summary-stops-length.tsv file.
  --prefix             Prefix for output files.
```

6) Combine alignment results, read mapping data, and incorporate functional annotations (emapper output)

Usage/help:

```bash
python3 tools/incorporate_annotations.py -h
usage: incorporate_annotations.py [-h] --annotations_file ANNOTATIONS_FILE --stops_length_file STOPS_LENGTH_FILE --disruptions_file DISRUPTIONS_FILE --coverage_file COVERAGE_FILE --variant_file VARIANT_FILE

Incorporate functional annotations and coverage data into the pseudogene disruptions summary file.

options:
  -h, --help            show this help message and exit
  --annotations_file ANNOTATIONS_FILE
                        Path to the functional annotations file (e.g., emapper output).
  --stops_length_file STOPS_LENGTH_FILE
                        Path to the <prefix>-summary-stops-length.tsv file.
  --disruptions_file DISRUPTIONS_FILE
                        Path to the pseudogene disruptions summary TSV file.
  --coverage_file COVERAGE_FILE
                        Path to the read mapping coverage file <prefix>-read-coverage.txt.
  --variant_file VARIANT_FILE
                        Path to the variant file <prefix>-filtered-variants.txt.
```

7) Run Pseudofinder (alternative/complimentary pseudogene prediction)

Example:

```bash
# Make a protein database for Pseudofinder (if using a local reference)
makeblastdb -in /path/to/reference/My_reference_NCBI_protein.faa -dbtype prot -parse_seqids

# Run Pseudofinder annotate
python3 pseudofinder.py annotate --genome prokka_out/My_MAG.gbk --outprefix My_MAG --database /path/to/reference/My_reference_NCBI_protein.faa --threads 16
```

8) Compare results of the manual pipeline with Pseudofinder output

Usage/help:

```bash
python3 tools/compare_manual_with_pseudofinder.py -h
usage: compare_manual_with_pseudofinder.py [-h] --final_summary_file FINAL_SUMMARY_FILE --pseudofinder_gff PSEUDOFINDER_GFF --annotations_file ANNOTATIONS_FILE --coverage_file COVERAGE_FILE --variant_file VARIANT_FILE

Compare your pseudogene results with Pseudofinder output.

options:
  -h, --help            show this help message and exit
  --final_summary_file FINAL_SUMMARY_FILE
                        Path to your results file.
  --pseudofinder_gff PSEUDOFINDER_GFF
                        Path to the Pseudofinder GFF file.
  --annotations_file ANNOTATIONS_FILE
                        Path to the functional annotations file (e.g., emapper output).
  --coverage_file COVERAGE_FILE
                        Path to the read mapping coverage file.
  --variant_file VARIANT_FILE
                        Path to the variant file.
```

---

## Notes & Requirements

- Ensure required tools are installed: Prokka, Bowtie2, samtools, bcftools, BLAST (makeblastdb/blastn), Python 3 with dependencies listed in `requirements.txt`.
- The pipeline expects Prokka output (`.ffn`, `.faa`, `.gbk`) and a reference amino acid FASTA when comparing lengths.
- Jupyter notebook in `metabolism/` contains downstream analysis and integration examples.

---
