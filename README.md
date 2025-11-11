Pipeline to detect candidate pseudogenes in MAGs.

Directory structure:
input/
input/My_reference_cds.fna
input/My_MAG_genomic.fa
input/prokka_out/
input/prokka_out/My_MAG.ffn
input/prokka_out/My_MAG.faa
input/prokka_out/emapper/
input/prokka_out/emapper/My_MAG.emapper.annotations
tools/
tools/detect_stops_compare_length.py
tools/remap_coverage_variants.sh
tools/compare_pseudogenes_with_reference.sh
tools/incorporate_annotations.py
tools/compare_manual_with_pseudofinder.py
Pseudofinder/
Pseudofinder/My_MAG_pseudos.gff
metabolism/
metabolism/pseudogenes_metabolism_integration.ipynb
Files:
requirements.txt: Python dependencies (Biopython, pandas, numpy)
pseudogenes_metabolism_workflow.txt: more detailed workflow
pseudogenes_metabolism_integration.ipynb: Jupyter Notebok for Data Analysis

Quick example:

1) Annotate with Prokka to get CDS and proteins:

prokka --outdir prokka_out --prefix My_MAG --addgenes --addmrna --locustag MYMAG --compliant --centre mRs --genus Genus --metagenome --evalue 1e-05 My_MAG_genomic.fa

2) Annotate with Emapper

emapper.py --cpu 32 \
          -m diamond  \
          -i My_MAG_proteins.faa  \
          -o My_MAG  \
          --output_dir ./emapper  \
          --report_orthologs --decorate_gff yes

3) Detect internal stops and compare lengths MAG vs reference proteins (provide a reference amino acid FASTA)

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

4)  Validate pseudogene candidates using read mapping and variant calls with remap_to_coverage_variants.sh
- Run remap of the MAG (Prokka output .ffn) agains the source metagenome using Bowtie2 and samtools
- Run variant calling with bcftools
- get read coverage and filter variants by quality

./tools/remap_coverage_variants.sh -h
Usage: remap_coverage_variants.sh --ffn_file <path> --metagenome_file <path> --prefix <prefix>

Arguments:
  --ffn_file           Path to the nucleotide CDS FASTA file (e.g., Prokka .ffn).
  --metagenome_file    Path to the metagenome file.
  --prefix             Prefix for output files.

5) Analyze Alignments between the pseudogene candidates and the reference genome using blastn (genomic sequences)
- Truncated proteins (alignment length < 70–80% of reference).
- Frameshifts or internal stop codons (Look for high mismatch or gapopen values).
- Percent identity (High identity suggests conserved regions).
- Parse the BLAST output to identify pseudogenes with disruptions (e.g., truncated alignments, frameshifts).

./tools/compare_pseudogenes_with_reference.sh -h
Usage: compare_pseudogenes_with_reference.sh --ffn_file <path> --reference_file <path> --stops_length_file <path> --prefix <prefix>

Arguments:
  --ffn_file           Path to the nucleotide CDS FASTA file (e.g., Prokka .ffn).
  --reference_file     Path to the reference genome file.
  --stops_length_file  Path to the mRsS25-Chr1-summary-stops-length.tsv file.
  --prefix             Prefix for output files.

6) Combine alignment results, read mapping data and incorporate functional annotations (emapper output)

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

7) Run Pseudofinder tool
# python3 pseudofinder.py annotate --genome GENOME.GBF --outprefix PREFIX --database /PATH/TO/NR/nr --threads 16

makeblastdb -in /path/to/reference/My_reference_NCBI_protein.faa -dbtype prot -parse_seqids
python3 pseudofinder.py annotate --genome prokka_out/My_MAG.gbk --outprefix My_MAG --database /path/to/reference/My_reference_NCBI_protein.faa --threads 16

8) Compare results of manual method with Pseudofinder

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


