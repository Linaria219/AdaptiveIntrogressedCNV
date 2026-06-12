#!/usr/bin/env python3
import os
import pandas as pd

def load_selection_file(chrom, selection_dir):
    file_path = os.path.join(selection_dir, f"{chrom}.pbs.txt.outlier_pvalue")
    if not os.path.exists(file_path):
        print(f"Warning: Selection file does not exist: {file_path}")
        return None
    
    try:
        data_start_line = 0
        expected_columns = 6
        
        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                parts = line.strip().split('\t')
                if len(parts) == expected_columns:
                    data_start_line = i
                    break
        
        if data_start_line == 0:
            df = pd.read_csv(file_path, sep='\t', header=0)
        else:
            df = pd.read_csv(file_path, sep='\t', skiprows=data_start_line, header=None)
            
            if df.shape[1] == expected_columns:
                df.columns = ['CHROM', 'POS', 'PBS', 'ZSCORE', 'PVALUE', 'STATUS']
            else:
                print(f"Warning: File {file_path} has unexpected number of columns: {df.shape[1]}")
                return None
        
        print(f"  Successfully loaded selection file {os.path.basename(file_path)} with {df.shape[0]} rows")
        return df
    except Exception as e:
        print(f"Error reading selection file {file_path}: {e}")
        return None

def load_introgression_file(chrom, introgression_dir, population):
    file_path = os.path.join(introgression_dir, f"{chrom}.msea.{population}.fd.noZscore.sis")
    if not os.path.exists(file_path):
        print(f"Warning: Introgression file does not exist: {file_path}")
        return None
    
    try:
        data_start_line = 0
        expected_columns = 6
        
        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                parts = line.strip().split('\t')
                if len(parts) == expected_columns:
                    data_start_line = i
                    break
        
        if data_start_line == 0:
            df = pd.read_csv(file_path, sep='\t', header=0)
        else:
            df = pd.read_csv(file_path, sep='\t', skiprows=data_start_line, header=None)
            
            if df.shape[1] == expected_columns:
                df.columns = ['CHROM', 'START', 'END', 'FD', 'PVALUE', 'STATUS']
            else:
                print(f"Warning: File {file_path} has unexpected number of columns: {df.shape[1]}")
                return None
        
        numeric_columns = ['START', 'END', 'FD', 'PVALUE']
        for col in numeric_columns:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        print(f"  Successfully loaded introgression file {os.path.basename(file_path)} with {df.shape[0]} rows")
        return df
    except Exception as e:
        print(f"Error reading introgression file {file_path}: {e}")
        return None

def find_max_pbs_in_region(selection_df, chrom, start, end):
    """Find the row with the maximum PBS value within a CNV region"""
    if selection_df is None or selection_df.empty:
        return "NA", "NA"
    
    try:
        if 'POS' not in selection_df.columns:
            if selection_df.shape[1] >= 2:
                selection_df = selection_df.copy()
                selection_df['POS'] = pd.to_numeric(selection_df.iloc[:, 1], errors='coerce')
        
        if 'PBS' not in selection_df.columns:
            if selection_df.shape[1] >= 3:
                selection_df = selection_df.copy()
                selection_df['PBS'] = pd.to_numeric(selection_df.iloc[:, 2], errors='coerce')
        
        in_region = selection_df[
            (selection_df['POS'] >= start) & 
            (selection_df['POS'] <= end) & 
            (pd.notna(selection_df['PBS']))
        ]
        
        if in_region.empty:
            return "NA", "NA"
        
        max_pbs_idx = in_region['PBS'].idxmax()
        max_pbs_row = selection_df.loc[max_pbs_idx]
        
        if selection_df.shape[1] >= 5:
            p_value = max_pbs_row.iloc[4] if 'PVALUE' not in max_pbs_row.index else max_pbs_row['PVALUE']
        else:
            p_value = "NA"
        
        if selection_df.shape[1] >= 6:
            a_status = max_pbs_row.iloc[5] if 'STATUS' not in max_pbs_row.index else max_pbs_row['STATUS']
        else:
            a_status = "NA"
        
        return p_value, a_status
    except Exception as e:
        print(f"  Error finding maximum PBS value: {e}")
        return "NA", "NA"

def find_max_fd_in_overlap_region(intro_df, chrom, start, end, buffer=100000):
    """Find the row with the maximum fD value within an extended CNV region"""
    if intro_df is None or intro_df.empty:
        return "NA", "NA"
    
    expanded_start = max(0, start - buffer)
    expanded_end = end + buffer
    
    try:
        if 'START' not in intro_df.columns:
            if intro_df.shape[1] >= 2:
                intro_df = intro_df.copy()
                intro_df['START'] = pd.to_numeric(intro_df.iloc[:, 1], errors='coerce')
        
        if 'END' not in intro_df.columns:
            if intro_df.shape[1] >= 3:
                intro_df = intro_df.copy()
                intro_df['END'] = pd.to_numeric(intro_df.iloc[:, 2], errors='coerce')
        
        if 'FD' not in intro_df.columns:
            if intro_df.shape[1] >= 4:
                intro_df = intro_df.copy()
                intro_df['FD'] = pd.to_numeric(intro_df.iloc[:, 3], errors='coerce')
        
        overlapping = intro_df[
            (intro_df['START'] <= expanded_end) & 
            (intro_df['END'] >= expanded_start) & 
            (pd.notna(intro_df['FD'])) &
            (intro_df['FD'] >= 0) &
            (intro_df['FD'] <= 1)
        ]
        
        if overlapping.empty:
            return "NA", "NA"
        
        # Find the row with the maximum fD value
        max_fd_idx = overlapping['FD'].idxmax()
        max_fd_row = intro_df.loc[max_fd_idx]
        
        # Extract p-value and I status
        if intro_df.shape[1] >= 5:
            p_value = max_fd_row.iloc[4] if 'PVALUE' not in max_fd_row.index else max_fd_row['PVALUE']
        else:
            p_value = "NA"
        
        if intro_df.shape[1] >= 6:
            i_status = max_fd_row.iloc[5] if 'STATUS' not in max_fd_row.index else max_fd_row['STATUS']
        else:
            i_status = "NA"
        
        return p_value, i_status
    except Exception as e:
        print(f"  Error finding maximum fD value: {e}")
        return "NA", "NA"

def main():
    target_file = "path/to/your/CNV/file"  
    selection_dir = "path/to/PBS/files"
    introgression_dir = "path/to/fD/files"
    output_file = "path/to/output"

    if not os.path.exists(target_file):
        print(f"Error: Target file does not exist: {target_file}")
        return
    
    try:
        print(f"Reading target file: {os.path.basename(target_file)}")
        target_df = pd.read_csv(target_file, sep='\t', header=None)
        
        if target_df.shape[1] != 8:
            print(f"Warning: Target file column count does not match expectation. Expected 8 columns, got {target_df.shape[1]} columns")
            print(f"Continuing processing, but results may be inaccurate")
        
        print(f"Successfully read target file, {target_df.shape[0]} rows and {target_df.shape[1]} columns")
        
        # Ensure the first three columns are of the correct data type
        target_df[0] = target_df[0].astype(str)  # Chromosome number
        target_df[1] = pd.to_numeric(target_df[1], errors='coerce')  # Start coordinate
        target_df[2] = pd.to_numeric(target_df[2], errors='coerce')  # End coordinate
        
    except Exception as e:
        print(f"Error reading target file: {e}")
        return
    
    # Extract chromosome list
    chromosomes = target_df[0].unique()
    print(f"Found {len(chromosomes)} chromosomes: {sorted(chromosomes)}")
    
    # Group by chromosome
    chromosome_groups = target_df.groupby(0)

    additional_info = []
    
    # Process CNVs for each chromosome
    for chrom, chrom_df in chromosome_groups:
        print(f"\nProcessing chromosome {chrom}...")
        
        # Load reference files
        print(f"  Loading selection information file...")
        selection_df = load_selection_file(chrom, selection_dir)
        
        print(f"  Loading denisovan information file...")
        deni_df = load_introgression_file(chrom, introgression_dir, "deni")
        
        print(f"  Loading neanderthal information file...")
        nean_df = load_introgression_file(chrom, introgression_dir, "nean")
        
        # Count the number of CNVs for this chromosome
        chrom_cnv_count = chrom_df.shape[0]
        print(f"  Chromosome {chrom} has {chrom_cnv_count} CNVs to process")
        
        # Process each row for this chromosome
        for idx, row in chrom_df.iterrows():
            start = row[1]
            end = row[2]
            
            # 1. Supplement selection information
            pbs_p, pbs_a = find_max_pbs_in_region(selection_df, chrom, start, end)
            
            # 2. Supplement denisovan information
            deni_p, deni_i = find_max_fd_in_overlap_region(deni_df, chrom, start, end)
            
            # 3. Supplement neanderthal information
            nean_p, nean_i = find_max_fd_in_overlap_region(nean_df, chrom, start, end)
            
            # Store supplementary information
            additional_info.append({
                'index': idx,
                'pbs_p': pbs_p if pbs_p is not None else "NA",
                'pbs_a': pbs_a if pbs_a is not None else "NA",
                'deni_p': deni_p if deni_p is not None else "NA",
                'deni_i': deni_i if deni_i is not None else "NA",
                'nean_p': nean_p if nean_p is not None else "NA",
                'nean_i': nean_i if nean_i is not None else "NA"
            })
        
        # Display processing progress
        if chrom_cnv_count > 0:
            found_pbs = sum(1 for info in additional_info[-chrom_cnv_count:] if info['pbs_p'] != "NA")
            found_deni = sum(1 for info in additional_info[-chrom_cnv_count:] if info['deni_p'] != "NA")
            found_nean = sum(1 for info in additional_info[-chrom_cnv_count:] if info['nean_p'] != "NA")
            
            print(f"  Chromosome {chrom} processing complete:")
            print(f"    Found PBS information: {found_pbs}/{chrom_cnv_count} ({found_pbs/chrom_cnv_count*100:.1f}%)")
            print(f"    Found denisovan information: {found_deni}/{chrom_cnv_count} ({found_deni/chrom_cnv_count*100:.1f}%)")
            print(f"    Found neanderthal information: {found_nean}/{chrom_cnv_count} ({found_nean/chrom_cnv_count*100:.1f}%)")
    
    additional_info.sort(key=lambda x: x['index'])
    
    output_lines = []
    
    for i, row in target_df.iterrows():
        info = next((item for item in additional_info if item['index'] == i), None)
        
        if info:
            # Generate new row: original 8 columns + supplementary 6 columns
            new_row = list(row.values) + [
                str(info['pbs_p']), str(info['pbs_a']),
                str(info['deni_p']), str(info['deni_i']),
                str(info['nean_p']), str(info['nean_i'])
            ]
        else:
            # If no supplementary information is found, fill with NA
            new_row = list(row.values) + ["NA"] * 6
        
        output_lines.append('\t'.join(map(str, new_row)))
    
    try:
        with open(output_file, 'w') as f:
            for line in output_lines:
                f.write(line + '\n')
        
        print(f"\nProcessing complete! Results saved to: {output_file}")
        
        # Calculate the number of columns in the new file
        if output_lines:
            new_columns = len(output_lines[0].split('\t'))
            print(f"Original columns: {target_df.shape[1]}, New file columns: {new_columns}")
        
        # Statistic missing value
        total_rows = len(additional_info)
        if total_rows > 0:
            missing_counts = {
                'pbs_p': sum(1 for info in additional_info if info['pbs_p'] == "NA"),
                'pbs_a': sum(1 for info in additional_info if info['pbs_a'] == "NA"),
                'deni_p': sum(1 for info in additional_info if info['deni_p'] == "NA"),
                'deni_i': sum(1 for info in additional_info if info['deni_i'] == "NA"),
                'nean_p': sum(1 for info in additional_info if info['nean_p'] == "NA"),
                'nean_i': sum(1 for info in additional_info if info['nean_i'] == "NA")
            }
            
            print("\nMissing value statistics:")
            for key, count in missing_counts.items():
                print(f"  {key}: {count} rows ({count/total_rows*100:.1f}%)")
                
    except Exception as e:
        print(f"Error writing output file: {e}")

if __name__ == "__main__":
    main()
