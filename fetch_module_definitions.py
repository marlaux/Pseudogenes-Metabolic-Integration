#!/usr/bin/env python3

import requests
import pandas as pd
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Fetch KEGG Module definitions for a list of module IDs."
    )
    parser.add_argument(
        "--module-list",
        required=True,
        help="Path to the input file containing module IDs, one per line."
    )
    args = parser.parse_args()

    # Read the input file containing module_ids
    module_ids_file = args.module_list
    try:
        with open(module_ids_file, 'r') as f:
            module_ids = [line.strip() for line in f.readlines()]
    except FileNotFoundError:
        print(f"Error: File '{module_ids_file}' not found.")
        sys.exit(1)

    # Initialize an empty list to store module data
    module_data = []

    # Loop through each module ID and fetch its definition
    for module_id in module_ids:
        url = f"https://rest.kegg.jp/get/md:{module_id}"
        response = requests.get(url)
        
        if response.status_code == 200:
            # Parse the module definition
            module_info = {}
            for line in response.text.splitlines():
                if line.startswith("ENTRY"):
                    module_info["ENTRY"] = line.split()[1]
                elif line.startswith("NAME"):
                    module_info["NAME"] = line[12:].strip()  # Extract the name
                elif line.startswith("CLASS"):
                    module_info["CLASS"] = line[12:].strip()  # Extract the class
                elif line.startswith("PATHWAY"):
                    # Append all pathways (if multiple) as a comma-separated string
                    if "PATHWAY" in module_info:
                        module_info["PATHWAY"] += f", {line[12:].strip()}"
                    else:
                        module_info["PATHWAY"] = line[12:].strip()
            
            # Add the parsed module info to the list
            module_data.append(module_info)
        else:
            print(f"Failed to fetch module {module_id}. Status code: {response.status_code}")

    # Convert the list of module data into a DataFrame
    module_df = pd.DataFrame(module_data)
    module_df.rename(columns={"ENTRY": "module", "NAME": "module_name", "CLASS": "module_class", "PATHWAY": "module_pathway"}, inplace=True)

    # Save the DataFrame to a CSV file
    # Generate the output filename by replacing '_list.txt' with 'definitions.txt'
    output_file = module_ids_file.replace('_list.txt', '_definitions.txt')
    module_df.to_csv(output_file, index=False)
    print(f"Data saved to {output_file}")

    # Print the DataFrame
    print(module_df.head())

if __name__ == "__main__":
    main()
