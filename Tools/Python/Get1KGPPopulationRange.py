#!/usr/bin/env python3
import gzip
import os
import glob
from collections import defaultdict

def extract_sample_ids_from_vcf(vcf_path):
    """
    Extracting a list of sample IDs from a VCF file
    Read only the comment lines at the beginning of the file until the #CHROM line is found.
    """
    sample_ids = []
    
    try:
        with gzip.open(vcf_path, 'rt') as vcf_file:
            for line in vcf_file:
                if line.startswith('#CHROM'):
                    parts = line.strip().split('\t')
                    
                    if len(parts) >= 10:
                        sample_ids = parts[9:]
                        print(f"{len(sample_ids)} sample IDs were extracted from the VCF file.")
                        break
                    else:
                        print(f"Warning: #CHROM line format is invalid, insufficient columns: {len(parts)}")
                        return []
        
        return sample_ids
    
    except Exception as e:
        print(f"Error occurred while reading the VCF file: {e}")
        return []

def load_population_samples(pop_dir):
    """
    Load a list of sample IDs for 26 different user groups.
    Return a dictionary: {User Group ID: List of Sample IDs}
    """
    population_samples = {}
    
    # Find all .txt files
    pop_files = glob.glob(os.path.join(pop_dir, "*.txt"))
    
    for pop_file in pop_files:
        # Extract population ID from the filename (remove .txt extension)
        pop_id = os.path.basename(pop_file).replace('.txt', '')
        
        try:
            with open(pop_file, 'r') as f:
                # Read sample IDs, removing whitespace
                samples = [line.strip() for line in f if line.strip()]
                population_samples[pop_id] = samples
            
            print(f"  Loading population {pop_id}: {len(samples)} samples")
            
        except Exception as e:
            print(f"Error occurred while reading population file {pop_file}: {e}")
    
    print(f"Total populations loaded: {len(population_samples)}")
    return population_samples

def find_population_column_ranges(vcf_samples, population_samples):
    """
    Find the column ranges for each population based on the VCF sample order and population sample lists.
    Return a dictionary: {Population ID: [(Start Column, End Column), ...]}
    """
    # Create a mapping from samples to column indices (0-based)
    sample_to_index = {sample: i for i, sample in enumerate(vcf_samples)}
    
    # Store the column indices for each population
    pop_column_indices = defaultdict(list)
    
    # Find the population for each sample
    sample_to_pop = {}
    for pop_id, pop_samples in population_samples.items():
        for sample in pop_samples:
            if sample in sample_to_index:
                sample_to_pop[sample] = pop_id
            else:
                # Can record samples not found, but not handling here
                pass
    
    # Group samples by population
    for sample, index in sample_to_index.items():
        if sample in sample_to_pop:
            pop_id = sample_to_pop[sample]
            pop_column_indices[pop_id].append(index)
    
    # Sort the column indices for each population
    pop_column_ranges = {}
    for pop_id, indices in pop_column_indices.items():
        if not indices:
            continue
            
        # Sort the indices
        indices.sort()
        
        # Find the continuous ranges
        ranges = []
        start = indices[0]
        end = indices[0]
        
        for i in range(1, len(indices)):
            if indices[i] == end + 1:
                # Continuous
                end = indices[i]
            else:
                # Discontinuous, record the current range
                ranges.append((start, end))
                start = indices[i]
                end = indices[i]
        
        # Add the last range
        ranges.append((start, end))
        
        # Convert to 1-based column numbers (offset by 9 positions)
        # Original VCF samples start from column 10, corresponding to index 0
        # We want to convert to 1-based: index 0 -> column 1
        shifted_ranges = []
        for start_idx, end_idx in ranges:
            shifted_start = start_idx + 1  # Start from 1
            shifted_end = end_idx + 1      # Start from 1
            shifted_ranges.append((shifted_start, shifted_end))
        
        pop_column_ranges[pop_id] = shifted_ranges
    
    return pop_column_ranges

def format_range_string(ranges):
    """
    Format the range list into a string
    Example: [(1, 5), (7, 10)] -> "1-5,7-10"
    """
    range_strings = []
    for start, end in ranges:
        if start == end:
            range_strings.append(str(start))
        else:
            range_strings.append(f"{start}-{end}")
    
    return ",".join(range_strings)

def main():
    vcf_path = "path/to/1KGP.vcf.gz"
    pop_dir = "path/to/population/samples/"
    output_dir = "path/to/output/"
    
    if not os.path.exists(vcf_path):
        print(f"Error: VCF file does not exist: {vcf_path}")
        return
    
    if not os.path.exists(pop_dir):
        print(f"Error: Population directory does not exist: {pop_dir}")
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    output_file = os.path.join(output_dir, "1KGPpoprange.txt")
    
    print("Starting processing...")
    print(f"VCF file: {vcf_path}")
    print(f"Population file directory: {pop_dir}")
    print(f"Output file: {output_file}")
    print("=" * 50)
    
    # Stage 1: Extract sample IDs from VCF
    print("\n1. Extracting sample IDs from VCF file...")
    vcf_samples = extract_sample_ids_from_vcf(vcf_path)
    
    if not vcf_samples:
        print("Error: Failed to extract sample IDs from VCF file")
        return
    
    # Stage 2: Load population sample lists
    print("\n2. Loading population sample lists...")
    population_samples = load_population_samples(pop_dir)
    
    if not population_samples:
        print("Error: Failed to load any population samples")
        return
    
    # Stage 3: Find column ranges for each population
    print("\n3. Finding column ranges for each population...")
    pop_ranges = find_population_column_ranges(vcf_samples, population_samples)
    
    # Statistics
    found_pops = len(pop_ranges)
    total_pops = len(population_samples)
    
    print(f"Found {found_pops}/{total_pops} populations with samples in the VCF file")
    
    # Stage 4: Write to output file
    print("\n4. Writing to output file...")
    try:
        with open(output_file, 'w') as f:
            # Write header
            f.write("#Population_ID\tColumn_Ranges(1-based)\n")
            
            # Sort by population ID and write
            for pop_id in sorted(pop_ranges.keys()):
                ranges = pop_ranges[pop_id]
                range_str = format_range_string(ranges)
                f.write(f"{pop_id}\t{range_str}\n")
        
        print(f"Results saved to: {output_file}")
        
        # Display summary information
        print("\n" + "=" * 50)
        print("Processing complete!")
        print(f"VCF total sample count: {len(vcf_samples)}")
        print(f"Processed population count: {found_pops}")
        
        # Display column ranges for each population
        print("\nPopulation column range summary:")
        for pop_id in sorted(pop_ranges.keys()):
            ranges = pop_ranges[pop_id]
            range_str = format_range_string(ranges)
            sample_count = sum((end - start + 1) for start, end in ranges)
            print(f"  {pop_id}: Column range {range_str} (Total {sample_count} samples)")
        
        # Check if any populations were not found
        if found_pops < total_pops:
            missing_pops = set(population_samples.keys()) - set(pop_ranges.keys())
            print(f"\nWarning: {len(missing_pops)} populations not found in the VCF file:")
            for pop_id in sorted(missing_pops):
                print(f"  {pop_id}")
        
    except Exception as e:
        print(f"Error writing to output file: {e}")
        return

if __name__ == "__main__":
    main()
