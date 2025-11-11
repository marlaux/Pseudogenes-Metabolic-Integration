#!/bin/bash

# Function to display help
function display_help() {
    echo "Usage: compare_pseudogenes_with_reference.sh --ffn_file <path> --reference_file <path> --stops_length_file <path> --prefix <prefix>"
    echo
    echo "Arguments:"
    echo "  --ffn_file           Path to the nucleotide CDS FASTA file (e.g., Prokka .ffn)."
    echo "  --reference_file     Path to the reference genome file."
    echo "  --stops_length_file  Path to the mRsS25-Chr1-summary-stops-length.tsv file."
    echo "  --prefix             Prefix for output files."
    echo
    exit 1
}

# Parse input arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --ffn_file)
            ffn_file="$2"
            shift 2
            ;;
        --reference_file)
            reference_file="$2"
            shift 2
            ;;
        --prefix)
            prefix="$2"
            shift 2
            ;;
        --stops_length_file)
            stops_length_file="$2"
            shift 2
            ;;
        -h|--help)
            display_help
            ;;
        *)
            echo "Unknown argument: $1"
            display_help
            ;;
    esac
done

# Check if required arguments are provided
if [[ -z "$ffn_file" || -z "$reference_file" || -z "$stops_length_file" || -z "$prefix" ]]; then
    echo "Error: Missing required arguments."
    display_help
fi

# Set variables from arguments
MAG_FASTA="$ffn_file"
REFERENCE_FASTA="$reference_file"
SUMMARY_STOPS="$stops_length_file"
PREFIX="$prefix"

# Define output files and database names
BLAST_DB="${REFERENCE_FASTA}"
BLAST_OUTPUT="output/pseudogenes-candidates-vs-ref-blastn.tsv"

# Define temporary files for identifiers and pseudogene candidates
PSEUDOGENES_CANDIDATES=$(mktemp)
IDENTIFIERS=$(mktemp)

# Ensure output directory exists
mkdir -p output

# Ensure temporary files are cleaned up on exit
trap "rm -f $IDENTIFIERS $PSEUDOGENES_CANDIDATES" EXIT

# Get pseudogenes identifiers
cat "$SUMMARY_STOPS" | grep "True" | cut -f 1 > "$IDENTIFIERS"

# Extract sequences matching the identifiers
seqkit grep -f "$IDENTIFIERS" "$MAG_FASTA" -o "$PSEUDOGENES_CANDIDATES"
# Prepare Reference Genome Data
makeblastdb -in "$REFERENCE_FASTA" -dbtype nucl -parse_seqids

# Align Pseudogene Candidates to the Reference
blastn -query "$PSEUDOGENES_CANDIDATES" -db "$BLAST_DB" -out "$BLAST_OUTPUT" -outfmt 6 -evalue 1e-5

# Add header to the BLAST output
sed -i '1i qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' "$BLAST_OUTPUT"

# Check if the BLAST output file was created and is not empty
if [[ -f "$BLAST_OUTPUT" && -s "$BLAST_OUTPUT" ]]; then
    echo "BLAST output file $BLAST_OUTPUT created successfully."
    rm "$PSEUDOGENES_CANDIDATES"
else
    echo "Error: BLAST output file $BLAST_OUTPUT was not created or is empty."
    exit 1
fi

# Parse and summarize blast output
python3 <<EOF
import csv

# Input and output files
blast_output = "output/pseudogenes-candidates-vs-ref-blastn.tsv"
summary_output = "output/pseudogenes_candidates_blast_disruptions_summary.tsv"

# Thresholds
alignment_length_threshold = 0.7  # 70% of reference length
percent_identity_threshold = 80.0  # High identity threshold
mismatch_threshold = 5  # Example threshold for mismatches
gapopen_threshold = 1  # Example threshold for gaps

# Parse BLAST output
with open(blast_output, "r") as infile, open(summary_output, "w") as outfile:
    reader = csv.reader(infile, delimiter="\t")
    writer = csv.writer(outfile, delimiter="\t")
    
    # Skip the header in the input file
    next(reader, None)
    
    # Write header for the summary file
    writer.writerow(["query", "Subject_ID", "Percent_Identity", "Alignment_Length",
                     "Mismatches", "Gap_Opens", "Disruption"])
    
    for row in reader:
        qseqid, sseqid, pident, length, mismatch, gapopen, *_ = row
        pident = float(pident)
        length = int(length)
        mismatch = int(mismatch)
        gapopen = int(gapopen)
        
        # Initialize disruption list for each row
        disruption = []
        
        # Replace hardcoded reference length with actual reference length
        # Check and calculate reference length for reverse sequences
        if int(row[7]) > int(row[8]):
            reference_length = int(row[7]) - int(row[8]) + 1
        else:
            reference_length = int(row[8]) - int(row[7]) + 1  # Calculate reference length from sstart and send
        if length < alignment_length_threshold * reference_length:
            disruption.append("Truncated")
        if mismatch > mismatch_threshold:
            disruption.append("Frameshift/Internal Stop")
        if gapopen > gapopen_threshold:
            disruption.append("Frameshift")
        if pident < percent_identity_threshold:
            disruption.append("Low Identity")
        
        # Write to summary if there are disruptions
        if disruption:
            writer.writerow([qseqid, sseqid, pident, length, mismatch, gapopen, "; ".join(disruption)])
EOF

# Check if the BLAST output file was created and is not empty
if [[ -f "output/pseudogenes_candidates_blast_disruptions_summary.tsv" && -s "output/pseudogenes_candidates_blast_disruptions_summary.tsv" ]]; then
    echo "Summarized output file output/pseudogenes_candidates_blast_disruptions_summary.tsv created successfully."
else
    echo "Error: Summarized output file output/pseudogenes_candidates_blast_disruptions_summary.tsv was not created or is empty."
    exit 1
fi
