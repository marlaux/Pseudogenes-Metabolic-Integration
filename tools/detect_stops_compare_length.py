#!/usr/bin/env python3
"""Detect internal stop codons in CDS FASTA, compare lengths vs reference using blastp, and summarize.

Inputs: 
 - nucleotide CDS FASTA (e.g., from Prokka .ffn)
 - mag_proteins.faa (FASTA)
 - ref_proteins.faa (FASTA)

Output: combined TSV marking likely pseudogenes by occurrence of internal stop or alignment length off limits.
"""
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import subprocess
import pandas as pd
import tempfile
import argparse

def main(args):
    cds_fasta = args.cds_fasta
    mag_faa = args.mag_faa
    ref_faa = args.ref_faa
    prefix = args.prefix

    # Ensure the output directory exists
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Create a temporary file for stops_out_tsv with columns: query, has_internal_stop (0/1), first_internal_pos (nt), aa_length
    stops_out_tsv = tempfile.NamedTemporaryFile(delete=False, suffix="-internal-stops.tsv", dir=output_dir).name

    with open(stops_out_tsv, 'w') as out:
        out.write('query\thas_internal_stop\tfirst_internal_pos_nt\tcds_length\n')
        for rec in SeqIO.parse(cds_fasta, 'fasta'):
            seq = rec.seq
            # Check if sequence length is divisible by 3
            if len(seq) % 3 != 0:
                print(f"Warning: Sequence {rec.id} length ({len(seq)}) is not a multiple of 3. Padding with 'N'.")
                seq = seq + Seq('N' * (3 - len(seq) % 3))
            # translate with standard table, stop_symbol='*'
            aa = seq.translate(table=11, to_stop=False)
            aa_len = len(aa.rstrip('*'))
            # find internal '*' not at last position
            internal_pos = None
            for i, a in enumerate(aa):
                if a == '*':
                    # if not the terminal stop (i != len(aa)-1) or there are multiple stops
                    if i != len(aa) - 1:
                        internal_pos = i * 3 + 1
                        break
            has_internal = 1 if internal_pos is not None else 0
            # Use the original sequence length as nucleotide length
            nt_len = len(seq)
            out.write(f"{rec.id}\t{has_internal}\t{internal_pos or ''}\t{nt_len}\n")

    # Check if the output file was created and is not empty
    if os.path.exists(stops_out_tsv) and os.path.getsize(stops_out_tsv) > 0:
        print(f"Internal stop codons succesfully identified ({sum(1 for _ in open(stops_out_tsv)) - 1} rows)")
    else:
        print(f"Error during internal stop codons detection.")

    # Create a temporary file for length_out_tsv with the columns: query, mag_len, ref_id, ref_len, length_ratio
    length_out_tsv = tempfile.NamedTemporaryFile(delete=False, suffix="-compare-length.tsv", dir=output_dir).name

    # helper: load lengths
    def load_lengths(fasta):
        d = {}
        for r in SeqIO.parse(fasta, 'fasta'):
            d[r.id] = len(r.seq)
        return d

    mag_len = load_lengths(mag_faa)
    ref_len = load_lengths(ref_faa)

    # Always create the BLAST database and run BLASTP
    try:
        subprocess.run(['makeblastdb', '-in', ref_faa, '-dbtype', 'prot'], check=True)
    except Exception as e:
        print('Failed to make BLAST DB; please install BLAST or supply precomputed m8 file')
        sys.exit(1)

    blastp_m8 = os.path.join(output_dir, 'pseudogenes-candidates-vs-ref-blastp.tsv')
    subprocess.run(['blastp', '-query', mag_faa, '-db', ref_faa, '-out', blastp_m8, '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'], check=True)

    # parse best hit per query (keep top bitscore)
    best = {}
    with open(blastp_m8) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            q, s, pident, aln_len, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = parts
            bitscore = float(bitscore)
            if q not in best or best[q][0] < bitscore:
                best[q] = (bitscore, s)

    with open(length_out_tsv, 'w') as out:
        out.write('query\tref_id\taa_mag_len\taa_ref_len\tlength_ratio\n')
        for q in mag_len:
            s = best[q][1] if q in best else ''
            magl = mag_len[q]
            refl = ref_len.get(s, '') if s else ''
            ratio = ''
            if refl:
                try:
                    ratio = float(magl) / float(refl)
                except Exception:
                    ratio = ''
            out.write(f"{q}\t{s}\t{magl}\t{refl}\t{ratio}\n")

    # Check if the output file was created and is not empty
    if os.path.exists(length_out_tsv) and os.path.getsize(length_out_tsv) > 0:
        print(f"Alignment length comparison succesfully completed ({sum(1 for _ in open(length_out_tsv)) - 1} rows)")
    else:
        print(f"Error during alignment length comparison.")

    # create summarized output
    summary_out_file = os.path.join(output_dir, f"{prefix}-summary-stops-length.tsv")

    # Read the temporary files into DataFrames
    stops_df = pd.read_csv(stops_out_tsv, sep='\t')
    length_df = pd.read_csv(length_out_tsv, sep='\t')
        
    # merge on gene id (mag_id == gene_id)
    # Perform the merge and ensure empty values resulting from 'outer' are NaN
    merged = pd.merge(stops_df, length_df, on='query', how='outer').fillna(pd.NA)

    # simple heuristic: pseudogene if internal_stop or length_ratio < 0.75
    merged['length_ratio'] = pd.to_numeric(merged['length_ratio'], errors='coerce')
    merged['is_pseudogene'] = ((merged['has_internal_stop'] == 1) | (merged['length_ratio'] < 0.75))

    # Save the summary output file with NaN explicitly written
    merged.to_csv(summary_out_file, sep='\t', index=False, na_rep='NaN')

    # Check if the summary output file was created and is not empty
    if os.path.exists(summary_out_file) and os.path.getsize(summary_out_file) > 0:
        print(f"Internal stops and alignment length summary succesfully created.")
        print(f"Wrote {summary_out_file} with {sum(1 for _ in open(summary_out_file)) - 1} rows")

        # Remove temporary files
        os.remove(stops_out_tsv)
        os.remove(length_out_tsv)
    else:
        print(f"Error during summary creation.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Detect internal stop codons, compare lengths vs reference, and summarize likely pseudogenes."
    )
    parser.add_argument(
        "--cds_fasta",
        required=True,
        help="Path to the nucleotide CDS FASTA file (e.g., Prokka .ffn)."
    )
    parser.add_argument(
        "--mag_faa",
        required=True,
        help="Path to the MAG proteins FASTA file."
    )
    parser.add_argument(
        "--ref_faa",
        required=True,
        help="Path to the reference proteins FASTA file."
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for output files."
    )
  
    args = parser.parse_args()

    main(args)