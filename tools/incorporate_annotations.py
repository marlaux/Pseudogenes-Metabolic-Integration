import argparse
import pandas as pd
import os
import numpy as np

def incorporate_functional_annotations(annotations_file, stops_length_file, disruptions_file, coverage_file, variant_file):
    """
    Incorporate functional annotations and coverage data into the pseudogene disruptions summary file.

    Parameters:
        annotations_file (str): Path to the functional annotations file (e.g., emapper output).
        stops_length_file (str): Path to the mRsS25-Chr1-summary-stops-length.tsv file.
        disruptions_file (str): Path to the pseudogene blast disruptions summary TSV file.
        coverage_file (str): Path to the read mapping coverage file.
        variant_file (str): Path to the variant file.
    """
    # Define the output file path
    output_file = "output/final_pseudogene_candidates_disruptions_with_annotations.tsv"

    # Load the emapper annotations file, skipping the first 4 metadata lines
    annotations_df = pd.read_csv(
        annotations_file, 
        sep='\t', 
        skiprows=4, 
        usecols=['#query', 'Preferred_name','COG_category', 'EC', 'KEGG_ko', 'KEGG_Module','Description']
    )
    
    # Rename the '#query' column to 'query' for merging
    annotations_df.rename(columns={'#query': 'query'}, inplace=True)

    # Load the pseudogene disruptions summary file
    disruptions_df = pd.read_csv(disruptions_file, sep='\t')

    # Load the mRsS25-Chr1-summary-stops-length.tsv file
    stops_length_df = pd.read_csv(stops_length_file, sep='\t')

    # Filter rows where '[is_pseudogene] == True'
    stops_length_df = stops_length_df[stops_length_df['is_pseudogene'] == True]

    # Identify queries absent from pseudogenes_candidates_disruptions_summary.tsv
    absent_queries = stops_length_df[~stops_length_df['query'].isin(disruptions_df['query'])].copy()

    # Create the '[Significant_Disruption]' column based on conditions
    def determine_significant_disruption(row):
        if row['has_internal_stop'] == 1 and row['length_ratio'] < 0.7:
            return 'Premature Stop/Len out range'
        elif row['has_internal_stop'] == 1:
            return 'Premature Stop'
        elif row['length_ratio'] < 0.7:
            return 'Align length out range'
        return None

    # Use .loc to avoid SettingWithCopyWarning
    absent_queries.loc[:, 'Disruption'] = absent_queries.apply(determine_significant_disruption, axis=1)

    # Merge the processed absent queries with the disruptions DataFrame using an outer join
    disruptions_df = pd.merge(disruptions_df, absent_queries[['query', 'Disruption']], on=['query', 'Disruption'], how='outer').sort_values(by='query').drop_duplicates()

    # Ensure all empty values are explicitly set to NaN
    disruptions_df = disruptions_df.replace({None: np.nan})

    # Load the read mapping coverage file with three columns
    coverage_df = pd.read_csv(coverage_file, sep='\t', header=None, names=['query', 'Position', 'Coverage'])

    # Aggregate coverage values by Query_ID (e.g., calculate mean coverage)
    # This step ensures coverage is calculated only once and used directly
    coverage_df = coverage_df.groupby('query', as_index=False)['Coverage'].mean()

    # Ensure query columns in both DataFrames are strings
    #disruptions_df['query'] = disruptions_df['query'].astype(str)
    #coverage_df['query'] = coverage_df['query'].astype(str)

    # Calculate mean and standard deviation of coverage for intact genes (not in pseudogene candidates)
    intact_genes_df = coverage_df[~coverage_df['query'].isin(disruptions_df['query'])]
    intact_mean_coverage = intact_genes_df['Coverage'].mean()
    intact_std_coverage = intact_genes_df['Coverage'].std()

    # Merge the coverage data with the disruptions summary first
    disruptions_with_coverage_df = disruptions_df.merge(coverage_df, on='query', how='left')

    # Ensure all empty values are explicitly set to NaN
    disruptions_with_coverage_df = disruptions_with_coverage_df.replace({None: np.nan})

    # Add a new column 'Cov_comp' based on the coverage comparison
    disruptions_with_coverage_df['Cov_comp'] = disruptions_with_coverage_df['Coverage'].apply(
        lambda cov: 'Low Cov' if cov < (intact_mean_coverage - 1 * intact_std_coverage) else 'Mean Cov'
    )

    # Load the variants file to include DP and Variant_Type columns
    variant_columns = ['query', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
    variants_df = pd.read_csv(variant_file, sep='\t', names=variant_columns, header=None, dtype={'query': str})

    # Extract DP and Variant_Type from the INFO column using the correct logic from filter_variants.py
    def extract_dp(info):
        try:
            for term in info.split(';'):
                if term.startswith('DP='):
                    return int(term.split('=')[1])
        except (ValueError, AttributeError):
            return None
        return None

    def extract_variant_types(info):
        try:
            return info.split(';')[0]
        except AttributeError:
            return None

    # Apply the extraction functions to the variants DataFrame
    variants_df['DP'] = variants_df['INFO'].apply(extract_dp)
    variants_df['Variant_Type'] = variants_df['INFO'].apply(extract_variant_types)

    # Ensure Variant_Type contains valid strings and handle NaN values
    variants_df['Variant_Type'] = variants_df['Variant_Type'].fillna('None').astype(str)

    # Merge DP and Variant_Type columns into the disruptions_with_coverage_df
    disruptions_with_coverage_var_df = disruptions_with_coverage_df.merge(
        variants_df[['query', 'DP', 'Variant_Type']], on='query', how='left'
    )

    # Merge the annotations after the coverage data
    disruptions_with_coverage_var_anno_df = disruptions_with_coverage_var_df.merge(annotations_df, on='query', how='left').sort_values(by='query').drop_duplicates()

    # Ensure all empty values are explicitly set to NaN
    disruptions_with_coverage_var_anno_df = disruptions_with_coverage_var_anno_df.replace({None: np.nan})
    
    # Reorder columns to include Variant_Type
    disruptions_with_coverage_var_anno_df = disruptions_with_coverage_var_anno_df[['query', 'Disruption', 'Coverage', 'Cov_comp', 'Variant_Type', 'DP',
                     'Preferred_name', 'COG_category', 'EC', 'KEGG_ko', 'KEGG_Module', 'Description']]
    
    # Replace Null and empty values with a dash '-'
    disruptions_with_coverage_var_anno_df = disruptions_with_coverage_var_anno_df.fillna('-')
    
    # Ensure the Coverage column is numeric before rounding
    disruptions_with_coverage_var_anno_df['Coverage'] = pd.to_numeric(disruptions_with_coverage_var_anno_df['Coverage'], errors='coerce')
    disruptions_with_coverage_var_anno_df['Coverage'] = disruptions_with_coverage_var_anno_df['Coverage'].round(2)
    
    # Save the updated summary file
    disruptions_with_coverage_var_anno_df.to_csv(output_file, sep='\t', index=False)

    # Check if the output file is not empty before printing the success message
    if not disruptions_with_coverage_var_anno_df.empty:
        print(f"Functional annotations and coverage data incorporated and saved to {output_file}")
    else:
        print(f"Warning: The output file {output_file} is empty. Please check the input files and parameters.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Incorporate functional annotations and coverage data into the pseudogene disruptions summary file.")
    parser.add_argument("--annotations_file", required=True, help="Path to the functional annotations file (e.g., emapper output).")
    parser.add_argument("--stops_length_file", required=True, help="Path to the <prefix>-summary-stops-length.tsv file.")
    parser.add_argument("--disruptions_file", required=True, help="Path to the pseudogene disruptions summary TSV file.")
    parser.add_argument("--coverage_file", required=True, help="Path to the read mapping coverage file <prefix>-read-coverage.txt.")
    parser.add_argument("--variant_file", required=True, help="Path to the variant file <prefix>-filtered-variants.txt.")

    args = parser.parse_args()

    incorporate_functional_annotations(
        disruptions_file=args.disruptions_file,
        annotations_file=args.annotations_file,
        stops_length_file=args.stops_length_file,
        coverage_file=args.coverage_file,
        variant_file=args.variant_file
    )