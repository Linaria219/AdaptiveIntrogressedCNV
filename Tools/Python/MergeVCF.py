#!/usr/bin/env python3
"""
VCF merge script: Used to merge VCF files from two population groups.
"""

import gzip
import os
import sys
from typing import List, Dict, Tuple, Optional

def parse_arguments():
    if len(sys.argv) != 2:
        print("Usage: python script.py <chromosome>")
        print("Example: python script.py 1")
        sys.exit(1)
    
    chrom = sys.argv[1]
    if not chrom.isdigit() or not (1 <= int(chrom) <= 22):
        print("Error: Chromosome number must be a digit between 1 and 22")
        sys.exit(1)
    
    return chrom

def get_file_paths(chrom: str) -> Tuple[str, str, str]:
    # Population 1
    pop1_dir = "path/to/pop1/vcf/files"
    pop1_file = os.path.join(pop1_dir, f"chr{chrom}.phased.vcf.gz")
    
    # Population 2
    kgp_dir = "path/to/pop1/vcf/files"
    kgp_file = os.path.join(kgp_dir, f"chrom_{chrom}.vcf.gz")
    
    output_file = f"chr{chrom}.pop1.pop2.merge.vcf"
    
    if not os.path.exists(pop1_file):
        print(f"Error: pop1 file does not exist: {pop1_file}")
        sys.exit(1)
    if not os.path.exists(kgp_file):
        print(f"Error: pop2 file does not exist: {kgp_file}")
        sys.exit(1)
    
    return pop1_file, kgp_file, output_file

def read_vcf_header(filename: str) -> Tuple[List[str], List[str], int]:
    header_lines = []
    sample_ids = []
    data_start_line = 0
    
    open_func = gzip.open if filename.endswith('.gz') else open
    mode = 'rt' if filename.endswith('.gz') else 'r'
    
    with open_func(filename, mode) as f:
        for i, line in enumerate(f):
            if line.startswith('#'):
                header_lines.append(line)
                if line.startswith('#CHROM'):
                    # Parse sample ID
                    parts = line.strip().split('\t')
                    if len(parts) >= 10:
                        sample_ids = parts[9:]
                    data_start_line = i + 1
            else:
                break
                
    return header_lines, sample_ids, data_start_line

def is_snp_line(line: str) -> bool:
    """Check if the pop2 line is a SNP (contains VT=SNP)"""
    parts = line.strip().split('\t')
    if len(parts) >= 8:
        info_field = parts[7]
        return 'VT=SNP' in info_field
    return False

def process_chromosome(pop1_file: str, kgp_file: str, output_file: str):
    """Process the VCF file merging for a single chromosome"""
    print(f"Starting to process chromosome files...")
    print(f"pop1 file: {pop1_file}")
    print(f"pop2 file: {kgp_file}")
    print(f"Output file: {output_file}")
    
    print("Reading pop1 file header...")
    pop1_header, pop1_samples, pop1_data_start = read_vcf_header(pop1_file)
    
    print("Reading pop2 file header...")
    kgp_header, kgp_samples, kgp_data_start = read_vcf_header(kgp_file)
    
    print(f"pop1 samples: {len(pop1_samples)}")
    print(f"pop2 samples: {len(kgp_samples)}")
    
    output_header = []
    for line in pop1_header:
        if line.startswith('#CHROM'):
            # Merge sample IDs
            merged_samples = pop1_samples + kgp_samples
            chrom_line = line.strip().split('\t')
            # Keep the first 9 columns unchanged, append the merged sample IDs
            new_chrom_line = chrom_line[:9] + merged_samples
            output_header.append('\t'.join(new_chrom_line) + '\n')
        else:
            output_header.append(line)
    
    with open(output_file, 'w') as out_f:
        out_f.writelines(output_header)
    
    # Process data lines - use two-pointer algorithm to leverage sorting
    print("Starting to process data lines (using two-pointer algorithm)...")
    
    pop1_open = gzip.open if pop1_file.endswith('.gz') else open
    kgp_open = gzip.open if kgp_file.endswith('.gz') else open
    
    pop1_mode = 'rt' if pop1_file.endswith('.gz') else 'r'
    kgp_mode = 'rt' if kgp_file.endswith('.gz') else 'r'
    
    matched_count = 0
    processed_lines = 0
    
    with pop1_open(pop1_file, pop1_mode) as pop1_f, \
         kgp_open(kgp_file, kgp_mode) as kgp_f, \
         open(output_file, 'a') as out_f:
        
        # Skip header lines
        for _ in range(pop1_data_start):
            next(pop1_f)
        for _ in range(kgp_data_start):
            next(kgp_f)
        
        # Initialize two pointers
        pop1_line = pop1_f.readline()
        kgp_line = kgp_f.readline()
        
        # Two-pointer traversal
        while pop1_line and kgp_line:
            processed_lines += 1
            if processed_lines % 100000 == 0:
                print(f"Processed {processed_lines} lines, matched {matched_count} SNVs")
            
            pop1_parts = pop1_line.strip().split('\t')
            kgp_parts = kgp_line.strip().split('\t')
            
            if len(pop1_parts) < 2 or len(kgp_parts) < 2:
                # Skip lines with incorrect format
                if len(pop1_parts) < 2:
                    pop1_line = pop1_f.readline()
                if len(kgp_parts) < 2:
                    kgp_line = kgp_f.readline()
                continue
            
            try:
                pop1_pos = int(pop1_parts[1])
                kgp_pos = int(kgp_parts[1])
            except ValueError:
                # Position is not a number, skip
                if not pop1_parts[1].isdigit():
                    pop1_line = pop1_f.readline()
                if not kgp_parts[1].isdigit():
                    kgp_line = kgp_f.readline()
                continue
            
            # Compare positions
            if pop1_pos < kgp_pos:
                pop1_line = pop1_f.readline()
            elif pop1_pos > kgp_pos:
                # Skip non-SNP lines
                next_kgp_line = kgp_f.readline()
                while next_kgp_line:
                    kgp_parts_next = next_kgp_line.strip().split('\t')
                    if len(kgp_parts_next) >= 2:
                        try:
                            next_kgp_pos = int(kgp_parts_next[1])
                            if next_kgp_pos >= pop1_pos:
                                kgp_line = next_kgp_line
                                break
                        except ValueError:
                            pass
                    next_kgp_line = kgp_f.readline()
                else:
                    kgp_line = None
            else:  # Positions are equal
                # Check if the pop2 line is a SNP
                if is_snp_line(kgp_line):
                    # Merge lines: complete pop1 line + GT column from pop2 (starting from column 10)
                    merged_line = '\t'.join(pop1_parts) + '\t' + '\t'.join(kgp_parts[9:]) + '\n'
                    out_f.write(merged_line)
                    matched_count += 1
                
                # Move pointers
                pop1_line = pop1_f.readline()
                kgp_line = kgp_f.readline()
    
    print(f"Processing complete! Total matched SNVs: {matched_count}")

def main():
    chrom = parse_arguments()
    
    try:
        pop1_file, kgp_file, output_file = get_file_paths(chrom)
        process_chromosome(pop1_file, kgp_file, output_file)
        print(f"Successfully generated merged file: {output_file}")
    except Exception as e:
        print(f"Error occurred during processing: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
