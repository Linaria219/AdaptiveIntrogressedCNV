#!/usr/bin/env python3
"""
CNV Statistical Analysis Script
Function: Performs CNV statistical analysis and status determination on BED files
Input: BED file, column parameters, statistical mode
Output: Result file containing statistics and status
"""

import argparse
import sys
import math
from scipy import stats
import numpy as np
from collections import defaultdict

def parse_column_spec(column_spec):
    """Parse column specification string (e.g., '1,3,6-19') into a list of column indices"""
    columns = []
    if not column_spec:
        return columns
    
    for part in column_spec.split(','):
        if '-' in part:
            start, end = map(int, part.split('-'))
            columns.extend(range(start-1, end))  # transform to 0-based index
        else:
            columns.append(int(part)-1)  # transform to 0-based index
    
    return sorted(set(columns))  # remove duplicates and sort

def strict_round(value):
    """Strict rounding function (not rounding to even numbers)."""
    if value < 0:
        return -strict_round(-value)
    
    integer_part = int(value)
    fractional_part = value - integer_part
    
    if fractional_part < 0.5:
        return integer_part
    elif fractional_part > 0.5:
        return integer_part + 1
    else:
        # In the case of 0.5, check the integer part.
        return integer_part + 1 if integer_part % 2 == 0 else integer_part

def detect_number_format(values, sample_ratio=0.05):
    """Detect the format of numeric values (integers or floats)"""
    if not values:
        return 'float' 
    
    # Sampling and Detection
    sample_size = max(1, int(len(values) * sample_ratio))
    sample_indices = np.random.choice(len(values), min(sample_size, len(values)), replace=False)
    
    integer_count = 0
    total_count = 0
    
    for idx in sample_indices:
        try:
            val = float(values[idx])
            # Check if it is an integer: the value minus its rounded value equals 0
            if abs(val - round(val)) < 1e-10:
                integer_count += 1
            total_count += 1
        except (ValueError, TypeError):
            continue
    
    if total_count == 0:
        return 'float'
    
    # If the majority of the sample are integers, consider it as integer format
    return 'int' if integer_count / total_count > 0.9 else 'float'

def convert_format(values, to_format):
    """Convert values to the specified format"""
    if to_format == 'int':
        return [str(strict_round(float(x))) for x in values]
    else:
        return [str(float(x)) for x in values]

def log2_transform(values):
    """Transform values using log2, handling values less than or equal to 0"""
    transformed = []
    for val in values:
        try:
            num_val = float(val)
            if num_val <= 0:
                num_val = 1e-6  # ten to the negative sixth power
            transformed.append(math.log2(num_val))
        except (ValueError, TypeError):
            transformed.append(math.log2(1e-6))
    return transformed

def calculate_vst(target_values, ref_values):
    """Calculate the Vst statistic"""
    if not target_values or not ref_values:
        return 0.0
    
    # Log2 transformation
    target_log = log2_transform(target_values)
    ref_log = log2_transform(ref_values)
    all_log = target_log + ref_log
    
    if len(all_log) < 2:
        return 0.0
    
    # Calculate variance
    VT = np.var(all_log, ddof=1)  # Total variance
    V_target = np.var(target_log, ddof=1) if len(target_log) > 1 else 0
    V_ref = np.var(ref_log, ddof=1) if len(ref_log) > 1 else 0
    
    M, N = len(target_log), len(ref_log)
    VS = (M/(M+N))*V_target + (N/(M+N))*V_ref if (M+N) > 0 else 0
    
    # Avoid division by zero and negative values
    if VT == 0:
        return 0.0
    
    vst = (VT - VS) / VT
    return max(0.0, min(vst, 1.0))  # Limit to 0-1 range

def calculate_mwu_bonferroni(target_values, ref_values, total_rows):
    """Calculate the Bonferroni p-value for the MWU test"""
    if not target_values or not ref_values:
        return 1.0
    
    try:
        # Convert to floats
        target_float = [float(x) for x in target_values]
        ref_float = [float(x) for x in ref_values]
        
        # Mann-Whitney U test
        u_stat, p_value = stats.mannwhitneyu(target_float, ref_float, alternative='two-sided')
        
        # Bonferroni correction
        bonferroni_p = min(p_value * total_rows, 1.0)
        return bonferroni_p
        
    except Exception:
        return 1.0

def calculate_dmedian(target_values, ref_values):
    """Calculate the Dmedian statistic"""
    if not target_values or not ref_values:
        return 0.0
    
    try:
        # Round strictly to integers
        target_int = [strict_round(float(x)) for x in target_values]
        ref_int = [strict_round(float(x)) for x in ref_values]
        
        # Calculate median
        target_median = np.median(target_int) if target_int else 0
        ref_median = np.median(ref_int) if ref_int else 0
        
        return abs(float(target_median - ref_median))
        
    except Exception:
        return 0.0

def determine_status(vst, bonferroni_p, dmedian):
    """Determine status based on three statistics"""
    try:
        # Ensure numerical comparison is correct (handling scientific notation)
        vst_float = float(vst)
        bonferroni_float = float(bonferroni_p)
        dmedian_float = float(dmedian)
        
        if (vst_float > 0.1 and 
            bonferroni_float < 0.05 and 
            dmedian_float > 0.5):
            return "Stratified"
        else:
            return "Neutral"
    except (ValueError, TypeError):
        return "Neutral"

def process_bed_file(args):
    """Process the BED file"""
    # Parse column specifications
    info_cols = parse_column_spec(args.info)
    target_cols = parse_column_spec(args.target)
    ref_cols = parse_column_spec(args.ref)
    
    print(f"Info columns: {[x+1 for x in info_cols]}")
    print(f"Target columns: {[x+1 for x in target_cols]}")
    print(f"Reference columns: {[x+1 for x in ref_cols]}")
    
    # First scan: get total row count (for Bonferroni correction)
    print("Calculating total row count...")
    total_rows = 0
    with open(args.bed, 'r') as f:
        for line in f:
            if line.strip():
                total_rows += 1
    
    print(f"Total rows: {total_rows}")
    
    # Process the file
    results = []
    processed_rows = 0
    
    with open(args.bed, 'r') as f:
        for line in f:
            if not line.strip():
                continue
                
            processed_rows += 1
            if processed_rows % 1000 == 0:
                print(f"Processed {processed_rows} rows")
            
            parts = line.strip().split('\t')
            
            # Extract column data
            info_values = [parts[i] for i in info_cols if i < len(parts)]
            target_values = [parts[i] for i in target_cols if i < len(parts)]
            ref_values = [parts[i] for i in ref_cols if i < len(parts)]
            
            # Skip invalid rows
            if not target_values or not ref_values:
                continue
            
            # Detect and convert number formats (only in int mode)
            if args.mode == 'int':
                target_format = detect_number_format(target_values)
                ref_format = detect_number_format(ref_values)
                
                if target_format != ref_format:
                    if target_format == 'int' and ref_format == 'float':
                        ref_values = convert_format(ref_values, 'int')
                    elif target_format == 'float' and ref_format == 'int':
                        target_values = convert_format(target_values, 'int')
            
            # Calculate three statistics
            vst = calculate_vst(target_values, ref_values)
            bonferroni_p = calculate_mwu_bonferroni(target_values, ref_values, total_rows)
            dmedian = calculate_dmedian(target_values, ref_values)
            
            # Determine status
            status = determine_status(vst, bonferroni_p, dmedian)
            
            # Store results
            results.append({
                'info': info_values,
                'vst': vst,
                'bonferroni_p': bonferroni_p,
                'dmedian': dmedian,
                'status': status
            })
    
    return results

def write_output(results, output_file):
    """Write output file"""
    with open(output_file, 'w') as f:
        for result in results:
            # Info columns
            info_str = '\t'.join(str(x) for x in result['info'])
            # Statistics (formatted in scientific notation to ensure precision)
            vst_str = f"{result['vst']:.6e}"
            bonferroni_str = f"{result['bonferroni_p']:.6e}"
            dmedian_str = f"{result['dmedian']:.6e}"
            
            # Write row
            output_line = f"{info_str}\t{vst_str}\t{bonferroni_str}\t{dmedian_str}\t{result['status']}\n"
            f.write(output_line)

def main():
    parser = argparse.ArgumentParser(description='CNV Stratification Analysis Tool')
    parser.add_argument('--bed', required=True, help='Input BED file path')
    parser.add_argument('--info', required=True, help='Info column specifications (e.g., 1,3,6-19)')
    parser.add_argument('--target', required=True, help='Target column specifications (e.g., 1,3,6-19)')
    parser.add_argument('--ref', required=True, help='Reference column specifications (e.g., 1,3,6-19)')
    parser.add_argument('--mode', choices=['float', 'int'], default='float', 
                       help='Number processing mode (float: ignore format differences, int: detect and convert formats)')
    parser.add_argument('--output', required=True, help='Output file path')
    
    args = parser.parse_args()
    
    print("Starting to process BED file...")
    print(f"Mode: {args.mode}")
    
    try:
        # Process the file
        results = process_bed_file(args)
        
        # Write output
        write_output(results, args.output)
        
        # Statistic results
        stratified_count = sum(1 for r in results if r['status'] == 'Stratified')
        neutral_count = len(results) - stratified_count
        
        print(f"Processing complete! Total rows: {len(results)}")
        print(f"Stratified: {stratified_count} rows")
        print(f"Neutral: {neutral_count} rows")
        print(f"Results saved to: {args.output}")
        
    except Exception as e:
        print(f"Error occurred during processing: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
