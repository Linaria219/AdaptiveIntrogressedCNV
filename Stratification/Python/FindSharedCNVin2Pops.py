#!/usr/bin/env python3
"""
Function: Find shared CNV regions between target and reference populations and generate merged files
"""

import argparse
import sys
import os
from collections import defaultdict
import numpy as np

def validate_and_get_columns(file_path):
    """Validate column count consistency and return column count"""
    with open(file_path, 'r') as f:
        lines = [line.strip().split('\t') for line in f if line.strip()]
    
    if not lines:
        return 0
    
    # Check if all rows have the same number of columns
    first_line_cols = len(lines[0])
    for i, parts in enumerate(lines[1:], 2):
        if len(parts) != first_line_cols:
            raise ValueError(f"Error: File {file_path} line {i} has inconsistent column count. Expected {first_line_cols} columns, got {len(parts)} columns")
    
    return first_line_cols

def parse_float_value(value):
    """Convert value to float, handling integer and decimal formats"""
    try:
        return float(value)
    except ValueError:
        return None

def read_ref_file(ref_file, ref_col_count):
    """Read and filter reference file, unify number formats"""
    ref_data = []
    valid_line_count = 0
    
    with open(ref_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            parts = line.strip().split('\t')
            if len(parts) != ref_col_count:
                continue
            
            # Check if columns 5 to the last contain "."

            skip_line = False
            for i in range(4, len(parts)):
                if parts[i] == '.':
                    skip_line = True
                    break
            
            if skip_line:
                continue
            
            # Convert number formats and validate
            cn_values = []
            valid_cn = True
            for i in range(4, len(parts)):
                num_val = parse_float_value(parts[i])
                if num_val is None:
                    valid_cn = False
                    break
                # Unify format: convert integers to decimal format (e.g., 2 → 2.000)
                cn_values.append(f"{num_val:.3f}")
            
            if not valid_cn:
                continue
            
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            ref_data.append((chrom, start, end, cn_values))
            valid_line_count += 1
    
    print(f"Reference file filtering complete, valid lines: {valid_line_count}/{line_num}")
    return ref_data, valid_line_count

def read_target_file(target_file, target_col_count):
    """Read and filter target file, validate number formats"""
    target_data = []
    valid_line_count = 0
    
    with open(target_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            parts = line.strip().split('\t')
            if len(parts) != target_col_count:
                continue
            
            # Check if columns 4 to the last are valid numbers
            cn_values = []
            valid_cn = True
            for i in range(3, len(parts)):
                num_val = parse_float_value(parts[i])
                if num_val is None:
                    valid_cn = False
                    break
                cn_values.append(parts[i])  # Keep original format
            
            if not valid_cn:
                continue
            
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            target_data.append((chrom, start, end, cn_values))
            valid_line_count += 1
    
    print(f"target file filtering complete, valid lines: {valid_line_count}/{line_num}")
    return target_data, valid_line_count

def calculate_overlap(start1, end1, start2, end2):
    """Calculate the overlap ratio between two regions"""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    if overlap_end <= overlap_start:
        return 0.0

    overlap_length = overlap_end - overlap_start
    total_length = max(end1, end2) - min(start1, start2)

    return overlap_length / total_length if total_length > 0 else 0.0

def calculate_position_difference(start1, end1, start2, end2):
    """Calculate the position difference between two regions"""
    start_diff = abs(start1 - start2)
    end_diff = abs(end1 - end2)
    return max(start_diff, end_diff)

def count_non2_columns(cn_values):
    """Count the number of columns with CN values not equal to 2"""
    count = 0
    for cn_val in cn_values:
        try:
            if abs(float(cn_val) - 2.0) > 1e-6:  # floating-point comparison tolerance
                count += 1
        except ValueError:
            continue
    return count

def find_best_match(ref_matches, target_cn_values):
    """Pick the best matching reference row from multiple matches"""
    if not ref_matches:
        return None
    
    # Pick the row with the most columns having CN values not equal to 2
    best_match = None
    max_non2_count = -1
    
    for ref_start, ref_end, ref_cn_values in ref_matches:
        non2_count = count_non2_columns(ref_cn_values)
        
        if non2_count > max_non2_count:
            max_non2_count = non2_count
            best_match = (ref_start, ref_end, ref_cn_values)
        elif non2_count == max_non2_count and best_match:
            # If the number of non-2 columns is the same, choose the row with CN values closer to the target (adjust as needed)
            current_diff = sum(abs(float(target_cn_values[i]) - float(best_match[2][i])) 
                             for i in range(len(target_cn_values)))
            new_diff = sum(abs(float(target_cn_values[i]) - float(ref_cn_values[i])) 
                          for i in range(len(target_cn_values)))
            if new_diff < current_diff:
                best_match = (ref_start, ref_end, ref_cn_values)
    
    return best_match

def find_shared_cnvs(target_data, ref_data, RO, PD):
    """Find shared CNV regions, using intelligent RO adjustment strategy"""
    shared_regions = []
    
    # Group reference data by chromosome
    ref_by_chrom = defaultdict(list)
    for chrom, start, end, cn_values in ref_data:
        ref_by_chrom[chrom].append((start, end, cn_values))

    # Sort reference data for each chromosome
    for chrom in ref_by_chrom:
        ref_by_chrom[chrom].sort(key=lambda x: x[0])

    # Process each row in the target data
    processed_count = 0
    for target_chrom, target_start, target_end, target_cn_values in target_data:
        processed_count += 1
        if processed_count % 5000 == 0:
            print(f"Processed {processed_count} rows of target data")

        # Check if the current chromosome has reference data
        if target_chrom not in ref_by_chrom:
            shared_regions.append({
                'target_chrom': target_chrom,
                'target_start': target_start,
                'target_end': target_end,
                'target_cn_values': target_cn_values,
                'ref_cn_values': None,
                'status': 'UNIQ'
            })
            continue

        # Find matching CNVs in the reference data for the current chromosome
        found_shared = False
        shared_refs = []
        refs_for_chrom = ref_by_chrom[target_chrom]

        # Use binary search to quickly locate potentially overlapping reference CNVs
        left, right = 0, len(refs_for_chrom) - 1
        start_idx = 0

        while left <= right:
            mid = (left + right) // 2
            ref_start, ref_end, _ = refs_for_chrom[mid]

            if ref_end < target_start:
                left = mid + 1
            elif ref_start > target_end:
                right = mid - 1
            else:
                start_idx = mid
                while start_idx > 0 and refs_for_chrom[start_idx-1][1] >= target_start:
                    start_idx -= 1
                break

        # Check all potentially overlapping reference CNVs starting from start_idx
        for i in range(start_idx, len(refs_for_chrom)):
            ref_start, ref_end, ref_cn_values = refs_for_chrom[i]

            if ref_start > target_end:
                break

            # Calculate overlap ratio
            overlap_ratio = calculate_overlap(target_start, target_end, ref_start, ref_end)
            if overlap_ratio < RO:
                continue

            # Calculate position difference (if not disabled)
            if PD != -1:
                pos_diff = calculate_position_difference(target_start, target_end, ref_start, ref_end)
                if pos_diff > PD:
                    continue

            found_shared = True
            shared_refs.append((ref_start, ref_end, ref_cn_values))

        # Intelligent RO adjustment strategy
        current_ro = RO
        last_valid_matches = shared_refs.copy()
        
        while len(shared_refs) > 1 and current_ro + 0.05 <= 1.0:
            current_ro += 0.05
            new_matches = []
            
            for ref_start, ref_end, ref_cn_values in last_valid_matches:
                overlap_ratio = calculate_overlap(target_start, target_end, ref_start, ref_end)
                if overlap_ratio >= current_ro:
                    new_matches.append((ref_start, ref_end, ref_cn_values))
            
            if not new_matches:
                # No matches found after increasing RO, use the last valid matches
                shared_refs = last_valid_matches
                break
            elif len(new_matches) == 1:
                # Found a unique match
                shared_refs = new_matches
                break
            else:
                # Still multiple matches, continue to the next round
                last_valid_matches = new_matches
                shared_refs = new_matches

        # If still multiple matches, use selection strategy
        if len(shared_refs) > 1:
            best_match = find_best_match(shared_refs, target_cn_values)
            shared_refs = [best_match] if best_match else []
        
        # Record results
        if shared_refs:
            for ref_start, ref_end, ref_cn_values in shared_refs:
                shared_regions.append({
                    'target_chrom': target_chrom,
                    'target_start': target_start,
                    'target_end': target_end,
                    'target_cn_values': target_cn_values,
                    'ref_cn_values': ref_cn_values,
                    'status': 'SHARE'
                })
        else:
            shared_regions.append({
                'target_chrom': target_chrom,
                'target_start': target_start,
                'target_end': target_end,
                'target_cn_values': target_cn_values,
                'ref_cn_values': None,
                'status': 'UNIQ'
            })

    return shared_regions

def create_merged_file(shared_regions, output_file, ref_sample_count):
    """Create a merged CNV file"""
    with open(output_file, 'w') as f_out:
        for region in shared_regions:
            output_row = [
                region['target_chrom'],
                str(region['target_start']),
                str(region['target_end']),
                region['status']
            ]

            output_row.extend(region['target_cn_values'])

            if region['status'] == 'SHARE':
                output_row.extend(region['ref_cn_values'])
            else:
                output_row.extend(['2.000'] * ref_sample_count)

            f_out.write('\t'.join(output_row) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Find shared CNV regions between target and reference populations')
    parser.add_argument('--target-file', required=True, help='Path to the target sample CNV file')
    parser.add_argument('--ref-file', required=True, help='Path to the reference data file')
    parser.add_argument('--RO', type=float, default=0.6, help='Overlap ratio threshold')
    parser.add_argument('--PD', type=float, default=1000, help='Position difference threshold')

    args = parser.parse_args()

    # Check if files exist
    for file_path in [args.target_file, args.ref_file]:
        if not os.path.isfile(file_path):
            print(f"Error: File '{file_path}' does not exist")
            sys.exit(1)

    print("Starting file processing...")
    
    try:
        # Check column count consistency
        ref_col_count = validate_and_get_columns(args.ref_file)
        target_col_count = validate_and_get_columns(args.target_file)
        
        print(f"Reference file column count: {ref_col_count}, target file column count: {target_col_count}")
        
        if ref_col_count < 5:
            print("Error: Reference file requires at least 5 columns")
            sys.exit(1)
        if target_col_count < 4:
            print("Error: Target file requires at least 4 columns")
            sys.exit(1)
            
    except ValueError as e:
        print(e)
        sys.exit(1)

    # Read and filter files
    ref_data, ref_valid_count = read_ref_file(args.ref_file, ref_col_count)
    target_data, target_valid_count = read_target_file(args.target_file, target_col_count)
    
    if ref_valid_count == 0:
        print("Error: Reference file has no valid data rows")
        sys.exit(1)
    if target_valid_count == 0:
        print("Error: Target file has no valid data rows")
        sys.exit(1)

    print("Starting to find shared CNV regions...")
    print(f"RO threshold: {args.RO}, PD threshold: {args.PD}")

    # Find shared regions
    shared_regions = find_shared_cnvs(target_data, ref_data, args.RO, args.PD)

    share_count = len([r for r in shared_regions if r['status'] == 'SHARE'])
    unique_count = len([r for r in shared_regions if r['status'] == 'UNIQ'])

    print(f"Found {share_count} shared regions")
    print(f"Found {unique_count} unique regions")

    # Calculate reference sample count (from the first valid data row)
    ref_sample_count = len(ref_data[0][3]) if ref_data else 0

    # Create output file
    output_file = f"target.1889.1KGP.noEAS.kmean.RO{args.RO}.PD{args.PD}.unmask.bed"
    create_merged_file(shared_regions, output_file, ref_sample_count)

    print(f"Processing complete! Results saved to {output_file}")

if __name__ == "__main__":
    main()
