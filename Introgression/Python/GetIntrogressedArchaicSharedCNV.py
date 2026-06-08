#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CNV-fD Window Matching Script
Function: Based on the input CNV bed file, search for windows that meet the specified criteria in the corresponding chromosome's fD window file.
"""

import os
import sys
import argparse
import bisect
from typing import List, Tuple, Dict, Optional

def read_fd_file(fd_file_path: str) -> List[Tuple]:
    """
    Read the fD window file and retain only rows with fD values in the range [0,1].
    Adapt to a six-column format: chrom, start, end, fd_value, p_value, status

    Args:
        fd_file_path: Path to the fD window file

    Returns:
        List of valid fD windows, each element is (chrom, start, end, fd_value, p_value, status, line_data)
    """
    valid_windows = []

    try:
        with open(fd_file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                # Only need 6 columns: chrom, start, end, fd, p_value, status
                if len(fields) < 6:
                    continue

                chrom, start_str, end_str, fd_str, p_value_str, status = fields[0], fields[1], fields[2], fields[3], fields[4], fields[5]

                try:
                    start = int(start_str)
                    end = int(end_str)
                    fd_value = float(fd_str)
                    p_value = p_value_str  # Keep string format, which may contain scientific notation

                    # Only retain rows with fD values in the range [0,1]
                    if 0 <= fd_value <= 1:
                        valid_windows.append((chrom, start, end, fd_value, p_value, status, line.strip()))
                except ValueError:
                    continue

    except FileNotFoundError:
        print(f"Warning: fD file does not exist {fd_file_path}")
    except Exception as e:
        print(f"Error occurred while reading fD file {fd_file_path}: {str(e)}")

    return valid_windows

def find_overlapping_windows(cnv_start: int, cnv_end: int, fd_windows: List[Tuple],
                           search_range: int = 100000) -> List[Tuple]:
    """
    Find overlapping fD windows with the CNV or within a specified range

    Args:
        cnv_start: CNV start position
        cnv_end: CNV end position
        fd_windows: List of fD windows
        search_range: Search range (upstream and downstream extension)

    Returns:
        List of matching fD windows
    """
    overlapping = []

    # CNV extended search range
    search_start = max(0, cnv_start - search_range)
    search_end = cnv_end + search_range

    for window in fd_windows:
        w_start, w_end = window[1], window[2]

        # Check if the window overlaps with the extended CNV region
        if not (w_end < search_start or w_start > search_end):
            overlapping.append(window)

    return overlapping

def filter_fd_windows(fd_windows: List[Tuple], strict_mode: bool = False) -> List[Tuple]:
    """
    Filter fD windows based on specified criteria

    Args:
        fd_windows: List of fD windows
        strict_mode: Whether to use strict mode

    Returns:
        Filtered list of fD windows
    """
    filtered = []

    for window in fd_windows:
        fd_value, status = window[3], window[5]  # The 4th element is fd_value, the 6th is status

        if strict_mode:
            # Strict mode: status must be "introgressed"
            if status == "introgressed":
                filtered.append(window)
        else:
            # Non-strict mode: fd_value >= 0.2
            if fd_value >= 0.2:
                filtered.append(window)

    return filtered

def process_cnv_file(input_file: str, output_file: str, fd_base_path: str, strict_mode: bool = False):
    """
    Process CNV file and search for matching fD windows

    Args:
        input_file: Input CNV file path
        output_file: Output file path
        fd_base_path: fD file base path
        strict_mode: Whether to use strict mode
    """

    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    current_chrom = None
    fd_windows = []

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        outfile.write("#CHROM\tSTART\tEND\tINFO\tMATCHED_FD_WINDOWS\n")
        
        for line in infile:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue

            chrom, start_str, end_str, info = fields[0], fields[1], fields[2], fields[3]

            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                continue

            # If switching to a new chromosome, load the corresponding fD file
            if chrom != current_chrom:
                current_chrom = chrom
                # Construct file path based on chromosome number (remove possible 'chr' prefix)
                chrom_num = chrom.replace('chr', '')
                #fd_file_path = os.path.join(fd_base_path, f"chr{chrom_num}.msea.nean.fd.noZscore.sis")
                fd_file_path = os.path.join(fd_base_path, f"chr{chrom_num}.msea.deni.fd.noZscore.sis")
                fd_windows = read_fd_file(fd_file_path)
                print(f"Load loaded fD data for chromosome {chrom}: {len(fd_windows)} valid windows")

            # Find overlapping or within-range fD windows
            overlapping = find_overlapping_windows(start, end, fd_windows)

            # Filter windows based on mode
            if strict_mode:
                filtered = filter_fd_windows(overlapping, strict_mode=True)
            else:
                filtered = filter_fd_windows(overlapping, strict_mode=False)

            # Sort by fD value in descending order
            filtered.sort(key=lambda x: x[3], reverse=True)

            # Prepare output line
            output_fields = [chrom, start_str, end_str, info]

            if filtered:
                # Add matching fD window information (4 columns per window: fd_value, start, end, status)
                for window in filtered:
                    fd_value, w_start, w_end, status = window[3], window[1], window[2], window[5]
                    output_fields.extend([f"{fd_value:.6f}", str(w_start), str(w_end), status])
            else:
                # No matching windows, add 4 "."
                output_fields.extend([".", ".", ".", "."])

            # Write to output file
            outfile.write("\t".join(output_fields) + "\n")

    print(f"Processing complete! Results saved to: {output_file}")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='CNV-fD window matching tool')
    parser.add_argument('input_file', help='Input CNV bed file path')
    parser.add_argument('output_file', help='Output file path')
    parser.add_argument('--fd_path', default='path/to/fD/',
                       help='fD file base path (default: path/to/fD/)')
    parser.add_argument('--strict', action='store_true',
                       help='Use strict mode (6th column is "introgressed")')

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file does not exist {args.input_file}")
        sys.exit(1)

    # Check if fD file path exists
    if not os.path.exists(args.fd_path):
        print(f"Error: fD file path does not exist {args.fd_path}")
        sys.exit(1)

    print(f"Starting to process CNV file: {args.input_file}")
    print(f"fD file path: {args.fd_path}")
    print(f"Mode: {'Strict (introgressed)' if args.strict else 'Relaxed (fD≥0.2)'}")

    # Process CNV file
    process_cnv_file(args.input_file, args.output_file, args.fd_path, args.strict)

if __name__ == "__main__":
    main()
