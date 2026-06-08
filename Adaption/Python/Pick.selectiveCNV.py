#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CNV Region Selective Signal SNV Detection Script
Function: Detects regions containing selective signal SNVs within a CNV region based on a PBS p-value file.
Input: CNV bed file, PBS p-value file directory
Output: CNV region file containing selective SNVs
"""

import argparse
import os
import sys
from typing import List, Dict, Tuple, Optional

def parse_arguments():
    parser = argparse.ArgumentParser(description='Detecting selective signal SNV in CNV region')
    parser.add_argument('--bed', required=True, help='Input CNV bed file path')
    parser.add_argument('--pbs_dir', required=True, help='PBS p-value file directory path')
    parser.add_argument('--extendbp', type=int, default=0, help='BP extension for CNV region, default is 0')
    parser.add_argument('--minzscore', type=float, default=2.57, help='Zscore threshold, default is 2.57')
    
    return parser.parse_args()

def read_cnv_bed(bed_file: str) -> List[Tuple]:
    """
    Read CNV bed file

    Args:
        bed_file: bed file path

    Returns:
        List of CNV regions, each element is (chrom, start, end, tag)
    """
    cnv_regions = []
    try:
        with open(bed_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    print(f"Warning: Skipping line {line_num}, insufficient columns: {line.strip()}")
                    continue
                
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                tag = parts[3]
                
                cnv_regions.append((chrom, start, end, tag))
        
        print(f"Loaded {len(cnv_regions)} CNV regions from {bed_file}")
        return cnv_regions
        
    except Exception as e:
        print(f"Error reading CNV bed file: {e}")
        sys.exit(1)

def read_pbs_file(pbs_file: str) -> Dict[int, Tuple]:
    """
    Read PBS p-value file
    
    Args:
        pbs_file: PBS file path
        
    Returns:
        Dictionary mapping positions to (PBS, Zscore, p_value, tag)
    """
    snv_data = {}
    if not os.path.exists(pbs_file):
        return snv_data
    
    try:
        with open(pbs_file, 'r') as f:
            # Skip header line
            header = next(f)
            
            for line_num, line in enumerate(f, 2):
                if not line.strip():
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    print(f"Warning: Skipping line {line_num}, insufficient columns: {line.strip()}")
                    continue
                
                try:
                    pos = int(parts[1])
                    pbs = float(parts[2])
                    zscore = float(parts[3])
                    p_value = float(parts[4])
                    tag = parts[5]
                    
                    snv_data[pos] = (pbs, zscore, p_value, tag)
                except ValueError as e:
                    print(f"Warning: Skipping line {line_num}, data format error: {e}")
                    continue
        
        print(f"Loaded {len(snv_data)} SNV sites from {pbs_file}")
        return snv_data
        
    except Exception as e:
        print(f"Error reading PBS file {pbs_file}: {e}")
        return {}

def is_selective_snv(snv_data: Tuple, min_zscore: Optional[float]) -> bool:
    """
    Determine if SNV is a selective signal
    
    Args:
        snv_data: SNV data tuple (PBS, Zscore, p_value, tag)
        min_zscore: Zscore threshold, if None then use tag to determine
        
    Returns:
        Whether it is a selective signal
    """
    pbs, zscore, p_value, tag = snv_data
    
    if min_zscore is not None:
        return zscore > min_zscore
    else:
        return tag == "adaptive"

def find_snv_in_cnv_region(cnv_region: Tuple, snv_data: Dict[int, Tuple], extend_bp: int) -> List[Tuple]:
    """
    Find the selective signal SNV within the CNV region
    
    Args:
        cnv_region: CNV region (chrom, start, end, tag)
        snv_data: SNV data dictionary
        extend_bp: Extension base pairs
        
    Returns:
        List of selective SNVs, each element is (pos, pbs, zscore)
    """
    chrom, start, end, tag = cnv_region
    
    # Using extended region to find SNVs
    extended_start = max(0, start - extend_bp)
    extended_end = end + extend_bp
    
    selective_snvs = []
    
    for pos, data in snv_data.items():
        if extended_start <= pos <= extended_end:
            pbs, zscore, p_value, snv_tag = data
            selective_snvs.append((pos, pbs, zscore))
    
    # Sort by PBS values in descending order
    selective_snvs.sort(key=lambda x: x[1], reverse=True)
    
    return selective_snvs

def process_cnv_regions(cnv_regions: List[Tuple], pbs_dir: str, extend_bp: int, min_zscore: Optional[float]) -> List[Tuple]:
    """
    Process all CNV regions, find selective signal SNVs
    
    Args:
        cnv_regions: List of CNV regions
        pbs_dir: Directory containing PBS files
        extend_bp: Number of base pairs to extend the CNV region
        min_zscore: Minimum Z-score threshold
        
    Returns:
        Processed results list, each element is CNV region info + selective SNV list
    """
    results = []
    
    # Group CNV regions by chromosome
    chrom_cnvs = {}
    for cnv in cnv_regions:
        chrom = cnv[0]
        if chrom not in chrom_cnvs:
            chrom_cnvs[chrom] = []
        chrom_cnvs[chrom].append(cnv)
    
    # Process each chromosome
    for chrom, chrom_regions in chrom_cnvs.items():
        # Build PBS file path
        pbs_file = os.path.join(pbs_dir, f"{chrom}.pbs.txt.outlier_pvalue")
        
        if not os.path.exists(pbs_file):
            print(f"Warning: PBS file does not exist: {pbs_file}")
            continue
        
        # Read PBS data for the chromosome
        snv_data = read_pbs_file(pbs_file)
        
        if not snv_data:
            print(f"Warning: Chromosome {chrom} has no valid PBS data")
            continue
        
        # Process each CNV region for the chromosome
        for cnv_region in chrom_regions:
            selective_snvs = find_snv_in_cnv_region(cnv_region, snv_data, extend_bp)
            
            # Only keep regions that contain selective signal SNVs
            if selective_snvs:
                results.append((cnv_region, selective_snvs))
    
    return results

def write_output(results: List[Tuple], output_file: str):
    """
    Args:
        results: Processed results list
        output_file: Output file path
    """
    try:
        with open(output_file, 'w') as f:
            # Write header
            f.write("#CHROM\tSTART\tEND\tTAG\tSELECTIVE_SNV_INFO\n")
            
            for cnv_region, selective_snvs in results:
                chrom, start, end, tag = cnv_region
                
                # Write CNV region basic information
                f.write(f"{chrom}\t{start}\t{end}\t{tag}")
                
                # Write information for each selective SNV
                for pos, pbs, zscore in selective_snvs:
                    f.write(f"\t{pos}\t{pbs:.6f}\t{zscore:.6f}")
                
                f.write("\n")
        
        print(f"Results written to file: {output_file}")
        print(f"Total {len(results)} CNV regions containing selective signal SNVs found")
        
    except Exception as e:
        print(f"Error occurred while writing output file: {e}")
        sys.exit(1)

def main():
    args = parse_arguments()
    
    # Parameter verification
    if args.extendbp < 0:
        print("Error: Extend base pairs must be a non-negative integer")
        sys.exit(1)
    
    if args.minzscore is not None and args.minzscore <= 0:
        print("Error: Z-score threshold must be a positive number")
        sys.exit(1)
    
    print("Starting processing of CNV regions for selective signal detection")
    print(f"Input CNV file: {args.bed}")
    print(f"PBS file directory: {args.pbs_dir}")
    print(f"Extension base pairs: {args.extendbp}")
    print(f"Z-score threshold: {args.minzscore if args.minzscore is not None else 'Using adaptive tag'}")
    
    # Read CNV bed file
    cnv_regions = read_cnv_bed(args.bed)
    
    if not cnv_regions:
        print("Error: No valid CNV regions found")
        sys.exit(1)
    
    # Process CNV regions
    results = process_cnv_regions(cnv_regions, args.pbs_dir, args.extendbp, args.minzscore)
    
    if not results:
        print("No CNV regions containing selective signal SNVs found")
        return
    
    # Prepare output file name
    input_basename = os.path.basename(args.bed)
    zscore_str = "adaptive" if args.minzscore is None else f"{args.minzscore}"
    output_file = f"selective.ex{args.extendbp}.zs{zscore_str}.{input_basename}"
    
    # Write to output file
    write_output(results, output_file)
    
    # Statistics
    total_selective_snvs = sum(len(snvs) for _, snvs in results)
    avg_snvs_per_cnv = total_selective_snvs / len(results)
    
    print(f"Processing complete!")
    print(f"Statistics:")
    print(f"  - CNV regions containing selective signal SNVs: {len(results)}")
    print(f"  - Total selective signal SNVs: {total_selective_snvs}")
    print(f"  - Average number of SNVs per CNV region: {avg_snvs_per_cnv:.2f}")

if __name__ == "__main__":
    main()
