#!/usr/bin/env python3
"""
Thanks to Andres Abreu for this script (2025).
ROOT File Variable Extractor

This script extracts specific variables from ROOT files using uproot.
It can handle single ROOT files, directories containing multiple ROOT files, 
or .txt files containing a list of ROOT file paths (one per line).
Variables can be specified using wildcard patterns (* ? [abc]).
"""

import argparse
import os
import glob
import sys
import fnmatch
from pathlib import Path

try:
    import uproot
    import awkward as ak
except ImportError as e:
    print(f"Error: Required package not found. Please install uproot and awkward:")
    print("pip install uproot awkward")
    print("Optional: pip install tqdm  # for progress bars")
    sys.exit(1)

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False


def read_file_list(txt_file):
    """
    Read ROOT file paths from a .txt file
    
    Args:
        txt_file (str): Path to .txt file containing ROOT file paths (one per line)
        
    Returns:
        list: List of ROOT file paths
    """
    root_files = []
    with open(txt_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            if not os.path.exists(line):
                print(f"Warning: File not found (line {line_num}): {line}")
                continue
            
            if not line.endswith('.root'):
                print(f"Warning: Not a ROOT file (line {line_num}): {line}")
                continue
            
            root_files.append(line)
    
    if not root_files:
        raise ValueError(f"No valid ROOT files found in {txt_file}")
    
    return root_files


def find_root_files(input_path):
    """
    Find ROOT files from input path (file, directory, or .txt file list)
    
    Args:
        input_path (str): Path to ROOT file, directory containing ROOT files, or .txt file with ROOT file paths
        
    Returns:
        list: List of ROOT file paths
    """
    if os.path.isfile(input_path):
        if input_path.endswith('.root'):
            return [input_path]
        elif input_path.endswith('.txt'):
            return read_file_list(input_path)
        else:
            raise ValueError(f"Input file {input_path} must be a ROOT file (.root) or a text file (.txt)")
    elif os.path.isdir(input_path):
        root_files = glob.glob(os.path.join(input_path, "*.root"))
        if not root_files:
            raise ValueError(f"No ROOT files found in directory {input_path}")
        return sorted(root_files)
    else:
        raise ValueError(f"Input path {input_path} does not exist")


def expand_variable_patterns(variable_patterns, available_branches):
    """
    Expand wildcard patterns in variable names to match available branches
    
    Args:
        variable_patterns (list): List of variable names/patterns (may contain wildcards)
        available_branches (list): List of available branch names in the tree
        
    Returns:
        list: List of expanded variable names that match the patterns
    """
    expanded_variables = []
    
    for pattern in variable_patterns:
        if '*' in pattern or '?' in pattern or '[' in pattern:
            # This is a wildcard pattern
            matches = [branch for branch in available_branches if fnmatch.fnmatch(branch, pattern)]
            if matches:
                expanded_variables.extend(sorted(matches))
                print(f"Pattern '{pattern}' matched {len(matches)} branches: {matches[:5]}{'...' if len(matches) > 5 else ''}")
            else:
                print(f"Warning: Pattern '{pattern}' matched no branches")
        else:
            # This is a literal variable name
            if pattern in available_branches:
                expanded_variables.append(pattern)
            else:
                print(f"Warning: Variable '{pattern}' not found in available branches")
    
    # Remove duplicates while preserving order
    seen = set()
    unique_variables = []
    for var in expanded_variables:
        if var not in seen:
            seen.add(var)
            unique_variables.append(var)
    
    return unique_variables


def extract_variables(input_files, tree_name, variables, output_file):
    """
    Extract specified variables from ROOT files and save to output file
    
    Args:
        input_files (list): List of input ROOT file paths
        tree_name (str): Name of the tree in ROOT files
        variables (list): List of variable names to extract
        output_file (str): Output ROOT file path
    """
    print(f"Processing {len(input_files)} ROOT file(s)...")
    
    # Get the actual variables to extract (expand wildcards using first file)
    actual_variables = None
    
    # Use first file to determine the actual variables from patterns
    for file_path in input_files:
        try:
            with uproot.open(file_path) as file:
                if tree_name in file:
                    tree = file[tree_name]
                    available_branches = list(tree.keys())
                    actual_variables = expand_variable_patterns(variables, available_branches)
                    print(f"Expanded to {len(actual_variables)} variables: {actual_variables[:10]}{'...' if len(actual_variables) > 10 else ''}")
                    break
        except Exception as e:
            print(f"Warning: Could not read {file_path} for variable expansion: {str(e)}")
            continue
    
    if actual_variables is None:
        raise ValueError("Could not determine variables from any input file")
    
    if not actual_variables:
        raise ValueError("No variables matched the specified patterns")
    
    all_data = {}
    total_entries = 0
    
    # Initialize data structure
    for var in actual_variables:
        all_data[var] = []
    
    # Process each input file
    file_iterator = enumerate(input_files)
    if HAS_TQDM and len(input_files) > 1:
        file_iterator = tqdm(file_iterator, total=len(input_files), desc="Processing files", unit="file")
    
    for i, file_path in file_iterator:
        if not HAS_TQDM or len(input_files) <= 1:
            print(f"Processing file {i+1}/{len(input_files)}: {os.path.basename(file_path)}")
        
        try:
            with uproot.open(file_path) as file:
                # Check if tree exists
                if tree_name not in file:
                    print(f"Warning: Tree '{tree_name}' not found in {file_path}")
                    available_trees = [key for key in file.keys() if hasattr(file[key], 'num_entries')]
                    if available_trees:
                        print(f"Available trees: {available_trees}")
                    continue
                
                tree = file[tree_name]
                
                # Check if all variables exist in the tree
                available_branches = tree.keys()
                missing_vars = [var for var in actual_variables if var not in available_branches]
                
                if missing_vars:
                    print(f"Warning: Variables {missing_vars} not found in {file_path}")
                    print(f"Available branches: {list(available_branches)[:10]}...")  # Show first 10
                    continue
                
                # Read the specified variables
                data = tree.arrays(actual_variables, library="ak")
                entries_in_file = len(data[actual_variables[0]])
                total_entries += entries_in_file
                
                if not HAS_TQDM or len(input_files) <= 1:
                    print(f"  - Read {entries_in_file} entries")
                
                # Append data to accumulated arrays
                for var in actual_variables:
                    all_data[var].append(data[var])
                
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            continue
    
    if total_entries == 0:
        print("Error: No data was successfully read from any input files")
        return False
    
    # Concatenate all data
    print(f"\nConcatenating data from all files...")
    final_data = {}
    
    var_iterator = actual_variables
    if HAS_TQDM and len(actual_variables) > 10:
        var_iterator = tqdm(actual_variables, desc="Concatenating variables", unit="var")
    
    for var in var_iterator:
        if all_data[var]:  # Only concatenate if we have data
            final_data[var] = ak.concatenate(all_data[var])
        else:
            if not HAS_TQDM or len(actual_variables) <= 10:
                print(f"Warning: No data found for variable '{var}'")
    
    if not final_data:
        print("Error: No variables were successfully extracted")
        return False
    
    # Write to output file
    print(f"Writing {total_entries} entries to {output_file}...")
    
    try:
        with uproot.recreate(output_file) as output:
            output[tree_name] = final_data
        
        print(f"Successfully created {output_file}")
        print(f"Total entries: {total_entries}")
        print(f"Variables written: {len(final_data)} data variables")
        
        # Check for counter variables but leave file as-is since they're needed for jagged arrays
        with uproot.open(output_file) as check_file:
            if tree_name in check_file:
                all_branches = list(check_file[tree_name].keys())
                counter_branches = [b for b in all_branches if b.startswith('n') and b[1:] in final_data]
                
                if counter_branches:
                    print(f"Note: {len(counter_branches)} counter variables automatically created")
                    print("(Counter variables are required by ROOT format for variable-length arrays)")
        
        return True
        
    except Exception as e:
        print(f"Error writing output file: {str(e)}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Extract specific variables from ROOT files using uproot",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract variables from a single ROOT file
  python script.py -i input.root -t MyTree -v var1 var2 var3

  # Extract variables from all ROOT files in a directory
  python script.py -i /path/to/root/files/ -t Events -v pt eta phi -o filtered_events.root

  # Extract variables from ROOT files listed in a .txt file
  python script.py -i file_list.txt -t Events -v pt eta phi mass -o combined_output.root

  # Extract all electron and muon variables using wildcards
  python script.py -i data.root -t Events -v 'Electron_*' 'Muon_*' -o leptons.root

  # Extract variables matching specific patterns
  python script.py -i input.root -t tree -v 'jet_*' '*_mass' 'pt_[0-9]' -o selected.root

  # Specify custom output file
  python script.py -i data.root -t tree -v branch1 branch2 -o output_custom.root
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input ROOT file, directory containing ROOT files, or .txt file with ROOT file paths (one per line)'
    )
    
    parser.add_argument(
        '-t', '--tree',
        required=True,
        help='Name of the tree in the ROOT files'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='output.root',
        help='Output ROOT file name (default: output.root)'
    )
    
    parser.add_argument(
        '-v', '--variables',
        nargs='+',
        required=True,
        help='List of variable names/patterns to extract from the ROOT files. Supports wildcards: * ? [abc]'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    try:
        input_files = find_root_files(args.input)
        print(f"Found {len(input_files)} ROOT file(s) to process")
        
        # Check if output file already exists
        if os.path.exists(args.output):
            response = input(f"Output file '{args.output}' already exists. Overwrite? (y/n): ")
            if response.lower() not in ['y', 'yes']:
                print("Operation cancelled")
                return
        
        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(args.output)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        
        # Extract variables
        success = extract_variables(
            input_files=input_files,
            tree_name=args.tree,
            variables=args.variables,
            output_file=args.output
        )
        
        if success:
            print("\nVariable extraction completed successfully!")
        else:
            print("\nVariable extraction failed!")
            sys.exit(1)
            
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
