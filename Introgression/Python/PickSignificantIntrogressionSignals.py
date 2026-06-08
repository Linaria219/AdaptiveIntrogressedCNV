#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fD file p-value calculation script: Directly calculate p-values and mark significant introgression signals
"""

import sys
import numpy as np
import os
from scipy import stats

def main():
    if len(sys.argv) != 2:
        print("Usage: python fd_pvalue.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = input_file + ".noZscore.sis"

    data = []  
    valid_fd_data = [] 
    header_lines = []

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line.strip())
                continue
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            chrom, start, end, fd_str = fields[0], fields[1], fields[2], fields[3]
            try:
                fd = float(fd_str)
                # Only consider fd values in the range [0,1] for valid data and background statistics calculation
                if 0 <= fd <= 1:
                    valid_fd_data.append(fd)
                    data.append((chrom, start, end, fd, True))  # Mark as valid
                else:
                    data.append((chrom, start, end, fd, False))  # Mark as invalid
            except ValueError:
                # Lines that cannot be converted to float are also considered invalid
                data.append((chrom, start, end, fd_str, False))

    if not data:
        print("Error: No valid data found in the file.")
        sys.exit(1)

    # Calculate background statistics (only based on 0<=fd<=1 data)
    if valid_fd_data:
        fd_values = np.array(valid_fd_data)
        mean_fd = np.mean(fd_values)
        std_fd = np.std(fd_values, ddof=1) if len(fd_values) > 1 else 0
    else:
        mean_fd, std_fd = 0, 0

    # Prepare output file
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    results = []

    # Process each window
    for item in data:
        chrom, start, end, fd_value, is_valid = item

        if is_valid and 0 <= fd_value <= 1 and std_fd != 0:
            # Directly calculate p-values (based on normal distribution assumption, but not explicitly calculating Z-scores)
            # Calculate standardized scores for p-value calculation, but do not store Z-scores
            z_score = (fd_value - mean_fd) / std_fd
            # Use two-sided p-values
            p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))

            # Apply screening criteria (introgressed standard unchanged: corresponding p<0.05)
            if  p_value < 0.05:
                status = "introgressed"
            else:
                status = "neutral"

            # Format output: do not include Z-score column, retain p-value and status
            results.append(f"{chrom}\t{start}\t{end}\t{fd_value:.6f}\t{p_value:.6e}\t{status}\n")
        else:
            # For invalid data or zero standard deviation cases, fill p-value and status columns with "."
            results.append(f"{chrom}\t{start}\t{end}\t{fd_value}\t.\t.\n")

    # Write to output file
    with open(output_file, 'w') as out_file:
        # Write original header lines
        for header in header_lines:
            out_file.write(header + "\n")

        # If no header is present, add a new header (without Z_SCORE column)
        if not header_lines:
            out_file.write("#CHROM\tSTART\tEND\tFD\tP_VALUE\tSTATUS\n")

        # Write data lines
        out_file.writelines(results)

    print(f"Processing complete. Output file: {output_file}")
    print(f"Total windows: {len(data)}")
    print(f"Valid windows for background statistics (0<=FD<=1): {len(valid_fd_data)}")
    if valid_fd_data:
        print(f"Background FD statistics: Mean = {mean_fd:.6f}, Standard deviation = {std_fd:.6f}")

if __name__ == "__main__":
    main()
