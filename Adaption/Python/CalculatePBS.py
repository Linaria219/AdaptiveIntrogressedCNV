#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PBS calculation script: Find the shared SNVs among the three groups using an efficient three-pointer synchronized traversal algorithm to optimize memory usage.
"""

import os
import sys
import logging
import math
from typing import Dict, List, Tuple, Optional, Generator

# logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def parse_fst_line(line: str) -> Optional[Tuple[int, int, float]]:
    """
    Parse the Fst file lines, handling possible variations in chromosome number format.

    Args:

    line: File line content

    Returns:

    (Chromosome number, position, Fst value) or None (if invalid)
    """
    parts = line.strip().split()
    if len(parts) < 3:
        return None
    
    try:
        # Process chromosome number formats (which may contain the "chr" prefix).
        chrom_str = parts[0]
        if chrom_str.startswith('chr'):
            chrom_str = chrom_str[3:]  # Remove the "chr" prefix
        chrom = int(chrom_str)
        
        pos = int(parts[1])
        fst = float(parts[2])
        
        # Check if Fst value is within valid range
        if 0 <= fst <= 1:
            return chrom, pos, fst
    except (ValueError, TypeError):
        pass
    
    return None

def calculate_pbs(fst1: float, fst2: float, fst3: float) -> float:
    """
    Calculate PBS value
    
    Args:
        fst1: A vs B Fst value
        fst2: A vs C Fst value
        fst3: B vs C Fst value

    Returns:
        PBS value
    """
    # Handling the case where Fst=1 (avoiding log(0))
    safe_fst = lambda x: 0.999999 if x >= 1.0 else x
    
    # Calculate T value
    T1 = -math.log(1 - safe_fst(fst1))
    T2 = -math.log(1 - safe_fst(fst2))
    T3 = -math.log(1 - safe_fst(fst3))
    
    # Calculate PBS value
    pbs = (T1 + T2 - T3) / 2
    return pbs

def read_fst_positions(file_path: str) -> Dict[int, float]:
    """
    Read all valid SNV positions and Fst values from the Fst file.

    Args:
        file_path: Path to the Fst file

    Returns:
        Dictionary mapping positions to Fst values
    """
    positions = {}
    try:
        with open(file_path, 'r') as f:
            # Skip header line
            for line in f:
                if line.startswith('CHROM') or line.startswith('chrom'):
                    break
            
            # Read data lines
            for line in f:
                data = parse_fst_line(line)
                if data:
                    chrom, pos, fst = data
                    positions[pos] = fst
        
        logger.debug(f"{len(positions)} valid SNVs were read from the file {file_path}.")
        return positions
    except Exception as e:
        logger.error(f"Error reading Fst file {file_path}: {str(e)}")
        return {}

def find_common_snvs_optimized(file1: str, file2: str, file3: str) -> Dict[int, Tuple[float, float, float]]:
    """
    Optimized Common SNV Finding Algorithm

    Using set operations to ensure all common sites are found

    Args:

    file1: Path to the first Fst file

    file2: Path to the second Fst file

    file3: Path to the third Fst file

    Returns:

    Dictionary of locations to the three Fst values
    """
    logger.info("Start searching for common SNV sites...")
    
    # Read all valid SNV locations from three files
    positions1 = read_fst_positions(file1)
    positions2 = read_fst_positions(file2) 
    positions3 = read_fst_positions(file3)
    
    # Find the common location of the three files.
    common_positions = set(positions1.keys()) & set(positions2.keys()) & set(positions3.keys())
    
    # Collect the Fst value of common SNV
    valid_snvs = {}
    for pos in common_positions:
        fst1 = positions1[pos]
        fst2 = positions2[pos] 
        fst3 = positions3[pos]
        valid_snvs[pos] = (fst1, fst2, fst3)
    
    logger.info(f"Found {len(valid_snvs)} common valid SNV sites")
    return valid_snvs

def process_chromosome_fst_files(chrom: int, fst_dir: str, output_dir: str) -> None:
    """
    Process Fst files for a single chromosome and calculate PBS values

    Args:

    chrom: Chromosome number

    fst_dir: Directory containing Fst files

    output_dir: Output directory
    """
    logger.info(f"Start processing chromosome {chrom}")
    
    # Build file paths
    file1 = os.path.join(fst_dir, f"chr{chrom}.A.vs.B.weir.fst")
    file2 = os.path.join(fst_dir, f"chr{chrom}.A.vs.C.weir.fst") 
    file3 = os.path.join(fst_dir, f"chr{chrom}.B.vs.C.weir.fst")
    
    # Check if files exist
    for f in [file1, file2, file3]:
        if not os.path.exists(f):
            logger.error(f"File does not exist: {f}")
            return
    
    # Find common SNVs
    valid_snvs = find_common_snvs_optimized(file1, file2, file3)
    
    if not valid_snvs:
        logger.warning(f"Chromosome {chrom} has no common valid SNVs")
        return
    
    # Prepare output file
    output_file = os.path.join(output_dir, f"chr{chrom}.pbs.txt")
    with open(output_file, 'w') as out_f:
        # Write header
        out_f.write("#CHROM\tPOS\tPBS\n")
        
        # Calculate and write PBS values
        for pos, (fst1, fst2, fst3) in valid_snvs.items():
            pbs = calculate_pbs(fst1, fst2, fst3)
            # The output format is: chr + number
            out_f.write(f"chr{chrom}\t{pos}\t{pbs:.6f}\n")
    
    logger.info(f"Complete chromosome {chrom}: Write {len(valid_snvs)} PBS values to {output_file}")

def main():
    fst_dir = "path/to/fst/files"  # Update this path to the directory containing Fst files
    output_dir = "path/to/output/directory"  # Update this path to the desired output directory
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Creating output directory: {output_dir}")
    
    # Process all autosomes (1-22)
    chromosomes = range(1, 23)
    total_common_snvs = 0
    
    for chrom in chromosomes:
        process_chromosome_fst_files(chrom, fst_dir, output_dir)
        
        pbs_file = os.path.join(output_dir, f"chr{chrom}.pbs.txt")
        if os.path.exists(pbs_file):
            with open(pbs_file, 'r') as f:
                line_count = sum(1 for line in f) - 1  # Remove header line
            total_common_snvs += line_count
            logger.info(f"Chromosome {chrom} processing complete, common SNV count: {line_count}")
    
    logger.info(f"All chromosomes processed! Total common SNV sites found: {total_common_snvs}")

if __name__ == "__main__":
    main()
