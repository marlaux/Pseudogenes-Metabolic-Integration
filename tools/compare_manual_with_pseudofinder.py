import pandas as pd
import argparse

def compare_with_pseudofinder(final_summary_file, pseudofinder_gff, annotations_file, coverage_file, variant_file):
    """
    Compare your pseudogene results with Pseudofinder output.

    Parameters:
        final_summary_file (str): Path to your results file (e.g., final_pseudogene_disruptions_with_annotations.tsv).
        pseudofinder_gff (str): Path to the Pseudofinder GFF file (e.g., pseudo_pseudos.gff).
    """
    # Load your results
    final_summary = pd.read_csv(final_summary_file, sep='\t')

    # Load the emapper annotations file, skipping the first 4 metadata lines
    annotations_df = pd.read_csv(
        annotations_file, 
        sep='\t', 
        skiprows=4, 
        usecols=['#query', 'Preferred_name','COG_category', 'EC', 'KEGG_ko', 'KEGG_Module','Description']
    )
    # Rename the '#query' column to 'query' for merging
    annotations_df.rename(columns={'#query': 'query'}, inplace=True)

    # Parse the Pseudofinder GFF file
    # Explode lists of values in Old_Locus_Tag and process entries containing 'ign'
    pseudofinder_data = []
    with open(pseudofinder_gff, 'r') as gff:
        for line in gff:
            if not line.startswith("#"):
                parts = line.strip().split('\t')
                attributes = {key: value for key, value in 
                              [attr.split('=') for attr in parts[8].split(';') if '=' in attr]}
                old_locus_tags = attributes.get('old_locus_tag', '').split(',')
                for tag in old_locus_tags:
                    if 'ign' in tag:
                        tag_parts = tag.split('_')
                        tag = f"{tag_parts[0]}_{tag_parts[-1]}"
                    pseudofinder_data.append({
                        'query': tag,
                        'Note': attributes.get('note', '')
                    })

    pseudofinder_df = pd.DataFrame(pseudofinder_data)

    # Compare results
    comparison = final_summary.merge(
        pseudofinder_df, 
        on='query', 
        how='outer', 
        indicator=True
    )

    # Add a column to indicate the source of identification
    comparison['Identification_Source'] = comparison['_merge'].map({
        'both': 'Match',
        'left_only': 'Manual Only',
        'right_only': 'Pseudofinder Only'
    })

    # Drop the '_merge' column as it's no longer needed
    comparison.drop(columns=['_merge'], inplace=True)

    # Sort the comparison DataFrame alphabetically by Query_ID (Old_Locus_Tag)
    comparison = comparison.sort_values(by='query').drop_duplicates()

    # replace 'Pseudogene candidate. Reason(s):' in 'Note' column with nothing (reduce the length)
    comparison['Note'] = comparison['Note'].str.replace('Pseudogene candidate. Reason(s):', '', regex=False)

    # rename and reorder columns to have 'query' and 'Identification_Source' first
    comparison.rename(columns={'Note': 'Pseudofinder_reason'}, inplace=True)

    comparison = comparison[['query', 'Identification_Source', 'Pseudofinder_reason','Disruption', 'Coverage', 'Cov_comp', 'Variant_Type', 'DP',
                     'Preferred_name', 'COG_category', 'EC', 'KEGG_ko', 'KEGG_Module', 'Description']]
    
    comparison1 = comparison[comparison['Identification_Source'] != 'Pseudofinder Only'].copy()
    
    # Convert Categorical columns to strings to avoid issues with fillna
    for col in comparison1.select_dtypes(['category']).columns:
        comparison1[col] = comparison1[col].astype(str)

    # Ensure '-' is added as a category for Categorical columns before filling NaN values
    for col in comparison1.select_dtypes(include=['category']).columns:
        comparison1[col] = comparison1[col].cat.add_categories('-')

    # Replace NaN values with '-' to maintain proper tabular structure
    comparison1 = comparison1.drop_duplicates().fillna('-')
    
    # Save the updated comparison report for 'Manual Only' and 'Match' entries with NaN replaced by '-'
    output_file1 = "output/comparison_manual_pseudofinder_report.tsv"

    # Check if the output file is not empty before printing the success message
    if not comparison1.empty:
        print(f"Manual and Match coverage and variants incorporated with functional annotations and saved to {output_file1}")
        comparison1.to_csv(output_file1, sep='\t', index=False)
    else:
        print(f"Warning: The output file {output_file1} is empty. Please check the input files and parameters.")

    # Filter for 'Pseudofinder Only' queries
    pseudofinder_only_df = comparison[comparison['Identification_Source'] == 'Pseudofinder Only'][['query','Identification_Source','Pseudofinder_reason']].copy()
    
    # Load the coverage file
    coverage_df = pd.read_csv(coverage_file, sep='\t', header=None, names=['query', 'Position', 'Coverage'])
    # Aggregate coverage values by query (e.g., calculate mean coverage)
    coverage_mean = coverage_df.groupby('query', as_index=False)['Coverage'].mean()

    # Calculate mean and standard deviation of coverage for intact genes (not in pseudogene candidates)
    intact_genes_df = coverage_df[~coverage_df['query'].isin(comparison['query'])]
    intact_mean_coverage = intact_genes_df['Coverage'].mean()
    intact_std_coverage = intact_genes_df['Coverage'].std()
    # Merge the coverage data with the disruptions summary first
    pseudofinder_with_coverage_df = pd.merge(pseudofinder_only_df,coverage_mean, on='query', how='left')
    # Add a new column 'Cov_comp' based on the coverage comparison
    pseudofinder_with_coverage_df['Cov_comp'] = pseudofinder_with_coverage_df['Coverage'].apply(
        lambda cov: 'Low Cov' if cov < (intact_mean_coverage - 1 * intact_std_coverage) else 'Mean Cov'
    )

    # Ensure the Coverage column is numeric and then round
    pseudofinder_with_coverage_df['Coverage'] = pd.to_numeric(pseudofinder_with_coverage_df['Coverage'], errors='coerce')
    pseudofinder_with_coverage_df['Coverage'] = pseudofinder_with_coverage_df['Coverage'].round(2)

    # Load the variants file to include DP and Variant_Type columns, ignoring the file's header
    variant_columns = ['query', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
    variants_df = pd.read_csv(variant_file, sep='\t', names=variant_columns, header=0, dtype={'query': str})

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
    pseudofinder_with_coverage_var_df = pseudofinder_with_coverage_df.merge(
        variants_df[['query', 'Variant_Type', 'DP']], on='query', how='left'
    )

    # Merge annotations
    pseudofinder_with_coverage_var_anno = pd.merge(pseudofinder_with_coverage_var_df,annotations_df, on='query', how='left').drop_duplicates()
    pseudofinder_with_coverage_var_anno = pseudofinder_with_coverage_var_anno[['query', 'Identification_Source', 'Pseudofinder_reason','Coverage', 'Cov_comp', 'Variant_Type', 'DP',
                     'Preferred_name', 'COG_category', 'EC', 'KEGG_ko', 'KEGG_Module', 'Description']]
    
    # Convert Categorical columns to strings to avoid issues with fillna
    for col in pseudofinder_with_coverage_var_anno.select_dtypes(['category']).columns:
        pseudofinder_with_coverage_var_anno[col] = pseudofinder_with_coverage_var_anno[col].astype(str)

    # Ensure '-' is added as a category for Categorical columns before filling NaN values
    for col in pseudofinder_with_coverage_var_anno.select_dtypes(include=['category']).columns:
        pseudofinder_with_coverage_var_anno[col] = pseudofinder_with_coverage_var_anno[col].cat.add_categories('-')

    # Replace NaN values with '-' to maintain proper tabular structure
    pseudofinder_with_coverage_var_anno = pseudofinder_with_coverage_var_anno.drop_duplicates().fillna('-')

    # Save the updated comparison report for 'Pseudofinder Only' entries with NaN replaced by '-'
    output_file2 = "output/comparison_pseudofinder_only_report.tsv"
    # Check if the output file is not empty before printing the success message
    if not pseudofinder_with_coverage_var_anno.empty:
        print(f"Pseudofinder only coverage and variants incorporated with functional annotations and saved to {output_file2}")
        pseudofinder_with_coverage_var_anno.to_csv(output_file2, sep='\t', index=False)
    else:
        print(f"Warning: The output file {output_file2} is empty. Please check the input files and parameters.")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare your pseudogene results with Pseudofinder output.")
    parser.add_argument("--final_summary_file", required=True, help="Path to your results file.")
    parser.add_argument("--pseudofinder_gff", required=True, help="Path to the Pseudofinder GFF file.")
    parser.add_argument("--annotations_file", required=True, help="Path to the functional annotations file (e.g., emapper output).")
    parser.add_argument("--coverage_file", required=True, help="Path to the read mapping coverage file.")
    parser.add_argument("--variant_file", required=True, help="Path to the variant file.")

    args = parser.parse_args()

    compare_with_pseudofinder(
        final_summary_file=args.final_summary_file,
        pseudofinder_gff=args.pseudofinder_gff,
        annotations_file=args.annotations_file,
        coverage_file=args.coverage_file,
        variant_file=args.variant_file
    )