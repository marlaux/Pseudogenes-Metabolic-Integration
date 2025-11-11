#!/bin/bash

# Function to display help
function display_help() {
    echo "Usage: remap_coverage_variants.sh --ffn_file <path> --metagenome_file <path> --prefix <prefix>"
    echo
    echo "Arguments:"
    echo "  --ffn_file           Path to the nucleotide CDS FASTA file (e.g., Prokka .ffn)."
    echo "  --metagenome_file    Path to the metagenome file."
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
        --metagenome_file)
            metagenome_file="$2"
            shift 2
            ;;
        --prefix)
            prefix="$2"
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
if [[ -z "$ffn_file" || -z "$metagenome_file" || -z "$prefix" ]]; then
    echo "Error: Missing required arguments."
    display_help
fi

# Set variables from arguments
PROKKA_FFN="$ffn_file"
METAGENOMES_DIR="$metagenome_file"
PREFIX="$prefix"
OUTPUT_DIR="output"
BOWTIE_INDEX="$OUTPUT_DIR/${PREFIX}-index"
SORTED_BAM="$OUTPUT_DIR/${PREFIX}-sorted.bam"
VCF_OUTPUT="$OUTPUT_DIR/${PREFIX}-variants.vcf.gz"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Check if the sorted BAM file already exists
if [[ -f "$SORTED_BAM" ]]; then
    echo "$SORTED_BAM already exists. Skipping alignment and proceeding to VCF generation."
else
    # Build bowtie2 index
    bowtie2-build "$PROKKA_FFN" "$BOWTIE_INDEX"

    # Find metagenome files
    READ1=$(find "$METAGENOMES_DIR" -name "*1*.fastq" | head -n 1)
    READ2=$(find "$METAGENOMES_DIR" -name "*2*.fastq" | head -n 1)

    if [[ -z "$READ1" || -z "$READ2" ]]; then
        echo "Error: Could not find paired metagenome files in $METAGENOMES_DIR"
        exit 1
    fi

    # Align reads and process in a single chain
    bowtie2 -p 32 -x "$BOWTIE_INDEX" -1 "$READ1" -2 "$READ2" |\
        samtools view -@ 32 -b - |\
        samtools sort -@ 32 -o "$SORTED_BAM"

    # Index the sorted BAM file
    samtools index -@ 32 "$SORTED_BAM"

    # Check if the sorted BAM file was created and is not empty
    if [[ -f "$SORTED_BAM" && -s "$SORTED_BAM" ]]; then
        echo "Sorted BAM file $SORTED_BAM created successfully."
    else
        echo "Error: Sorted BAM file $SORTED_BAM was not created or is empty."
        exit 1
    fi
fi

# Check if the VCF file already exists
if [[ -f "$VCF_OUTPUT" ]]; then
    echo "$VCF_OUTPUT already exists. Skipping mpileup and proceeding to filtering."
else
    # Generate VCF file directly
    bcftools mpileup -f "$PROKKA_FFN" -Ou "$SORTED_BAM" | \
        bcftools call -mv -Oz -o "$VCF_OUTPUT"

    # Check if the VCF output file was created and is not empty
    if [[ -f "$VCF_OUTPUT" && -s "$VCF_OUTPUT" ]]; then
        echo "VCF file $VCF_OUTPUT created successfully."
    else
        echo "Error: VCF file $VCF_OUTPUT was not created or is empty."
        exit 1
    fi
fi

# Generate filtered variants file with QUAL >= 15
bcftools view -i 'QUAL>=15' "$VCF_OUTPUT" | grep -v '^#' | awk '
BEGIN { OFS="\t" }
{
    ref = $4; alt = $5;
    if (length(ref) == 1 && length(alt) == 1) {
        if ((ref == "A" && alt == "G") || (ref == "G" && alt == "A") || (ref == "C" && alt == "T") || (ref == "T" && alt == "C")) {
            $8 = "Transition;" $8;
        } else {
            $8 = "Transversion;" $8;
        }
    }
    print $0;
}' > "$OUTPUT_DIR/${PREFIX}-filtered-variants.txt"

# Generate read coverage file
samtools depth -a "$SORTED_BAM" > "$OUTPUT_DIR/${PREFIX}-read-coverage.txt"

# Check if the read coverage file was created and is not empty
if [[ -f "$OUTPUT_DIR/${PREFIX}-read-coverage.txt" && -s "$OUTPUT_DIR/${PREFIX}-read-coverage.txt" ]]; then
    echo "Read coverage file $OUTPUT_DIR/${PREFIX}-read-coverage.txt created successfully."
else
    echo "Error: Read coverage file $OUTPUT_DIR/${PREFIX}-read-coverage.txt was not created or is empty."
    exit 1
fi

# Check if the filtered variants file was created and is not empty
if [[ -f "$OUTPUT_DIR/${PREFIX}-filtered-variants.txt" && -s "$OUTPUT_DIR/${PREFIX}-filtered-variants.txt" ]]; then
    echo "Filtered variants file $OUTPUT_DIR/${PREFIX}-filtered-variants.txt created successfully."
else
    echo "Error: Filtered variants file $OUTPUT_DIR/${PREFIX}-filtered-variants.txt was not created or is empty."
    exit 1
fi

# Add header to the filtered variants file
echo -e "query\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | cat - "$OUTPUT_DIR/${PREFIX}-filtered-variants.txt" > temp && mv temp "$OUTPUT_DIR/${PREFIX}-filtered-variants.txt"

# Check if all output files are non-empty before printing the completion message
if [[ -f "$VCF_OUTPUT" && -s "$VCF_OUTPUT" && \
      -f "$OUTPUT_DIR/${PREFIX}-filtered-variants.txt" && -s "$OUTPUT_DIR/${PREFIX}-filtered-variants.txt" && \
      -f "$OUTPUT_DIR/${PREFIX}-read-coverage.txt" && -s "$OUTPUT_DIR/${PREFIX}-read-coverage.txt" ]]; then
    echo "Remap and variant calling completed. Outputs are in $OUTPUT_DIR."
else
    echo "Error: One or more output files are missing or empty."
    exit 1
fi
