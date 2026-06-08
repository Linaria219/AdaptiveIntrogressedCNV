#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PBS Outlier Analysis and Adaptive Signal Identification Script
Function: Analyzes PBS files based on Z-score outlier detection, calculates the outlier p-value for each site, and identifies adaptive signals.
Input: Tab-delimited file containing PBS values.
Output: New file with added outlier p-values and status columns.
"""

import sys
import numpy as np
import scipy.stats as stats

def main():
    if len(sys.argv) < 2:
        print("Usage: python pbs_outlier_analysis.py input_file [z_threshold]")
        print("Default z_threshold is 3.0 (corresponding to p < 0.0027 for normal distribution)")
        sys.exit(1)
    
    input_file = sys.argv[1]
    z_threshold = 3.0  # The default Z-value threshold corresponds to a standard normal distribution with p ~ 0.0027.
    
    if len(sys.argv) >= 3:
        try:
            z_threshold = float(sys.argv[2])
        except ValueError:
            print("Error: z_threshold must be a numeric value")
            sys.exit(1)
    
    try:
        with open(input_file, 'r') as f:
            header = next(f).strip()
            data = []
            positions = []
            for line_num, line in enumerate(f, 2):  # Start from the second line
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    print(f"Warning: Skipping invalid line {line_num}: {line.strip()}")
                    continue
                try:
                    chr_num = parts[0]
                    pos = int(parts[1])
                    pbs = float(parts[2])
                    # Only keep sites with PBS >= 0 for subsequent analysis
                    if pbs >= 0:
                        data.append((chr_num, pos, pbs))
                        positions.append(pos)
                except ValueError as e:
                    print(f"Warning: Skipping line {line_num} due to conversion error: {e}")
                    continue
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
    
    if not data:
        print("Error: No valid data found in the input file")
        sys.exit(1)
    
    print(f"Loaded {len(data)} records with PBS >= 0 from {input_file}")
    
    # Sort by coordinate (second column)
    sorted_data = sorted(data, key=lambda x: x[1])
    print("Sorted data by coordinate")
    
    # Extract PBS values for outlier analysis
    pbs_values = np.array([row[2] for row in sorted_data])
    
    # Calculate Z-score (outlier indicator) for PBS values
    # Z-score = (X - μ) / σ, where μ is the mean and σ is the standard deviation
    pbs_mean = np.mean(pbs_values)
    pbs_std = np.std(pbs_values)
    
    # Avoid division by zero error
    if pbs_std == 0:
        print("Warning: All PBS values are identical, standard deviation is 0")
        z_scores = np.zeros_like(pbs_values)
    else:
        z_scores = (pbs_values - pbs_mean) / pbs_std
    
    print(f"PBS statistics: mean = {pbs_mean:.6f}, std = {pbs_std:.6f}")
    print(f"Z-score threshold: {z_threshold} (higher values indicate stronger outlier signals)")
    
    # Calculate the p-values corresponding to each Z-score (two-tailed test)
    # The p-value represents the probability of observing a similar or more extreme deviation under random conditions
    p_values = 2 * (1 - stats.norm.cdf(np.abs(z_scores)))
    
    # Determine adaptive signals (based on Z-score threshold)
    statuses = []
    adaptive_count = 0
    for z in z_scores:
        if z >= z_threshold:  # Focus only on positive outliers (high PBS values).
            statuses.append("adaptive")
            adaptive_count += 1
        else:
            statuses.append("neutral")
    
    print(f"Found {adaptive_count} adaptive signals (Z-score >= {z_threshold})")
    
    output_file = input_file + ".outlier_pvalue"
    
    try:
        with open(output_file, 'w') as f:
            f.write("#CHROM\tPOS\tPBS\tZ_score\tp_value\tstatus\n")
            
            for i, row in enumerate(sorted_data):
                chr_num, pos, pbs = row
                z_score = z_scores[i]
                p_val = p_values[i]
                status = statuses[i]
                f.write(f"{chr_num}\t{pos}\t{pbs:.6f}\t{z_score:.6f}\t{p_val:.6e}\t{status}\n")
        
        print(f"Results written to {output_file}")
        print(f"Total records processed: {len(sorted_data)}")
        print(f"Adaptive signals: {adaptive_count} ({adaptive_count/len(sorted_data)*100:.2f}%)")
        print(f"Top 5 strongest adaptive signals:")
        
        # Show the strongest 5 adaptive signals
        adaptive_indices = [i for i, z in enumerate(z_scores) if z >= z_threshold]
        if adaptive_indices:
            adaptive_with_z = [(i, z_scores[i], pbs_values[i]) for i in adaptive_indices]
            adaptive_sorted = sorted(adaptive_with_z, key=lambda x: x[1], reverse=True)[:5]
            
            for i, (idx, z, pbs_val) in enumerate(adaptive_sorted):
                chr_num, pos, _ = sorted_data[idx]
                print(f"  {i+1}. {chr_num}:{pos} - PBS={pbs_val:.4f}, Z={z:.4f}")
        
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
