#!/usr/bin/env python3
# gc_correct_cn.py

import pandas as pd
import numpy as np
import argparse
import os


def main(diploid_file, input_file, output_file):
    # Step 1: Read the diploid region file
    diploid_df = pd.read_csv(diploid_file, sep="\t")

    # Step 2: Read target window
    infer_df = pd.read_csv(input_file, sep="\t", header=None,
                           names=["chr", "start", "end", "GC_content", "mean_depth"])

    # Step 3: For each infer_df, find the K_loess of the corresponding bin for GC_content
    def lookup_k(gc_value):
        match = diploid_df[(diploid_df["bin_start"] <= gc_value) & (gc_value < diploid_df["bin_end"])]
        if not match.empty:
            return match.iloc[0]["K_loess"]
        else:
            return np.nan

    infer_df["K_loess"] = infer_df["GC_content"].apply(lookup_k)

    # Step 4: Calculate correction depth
    infer_df["corrected_depth"] = infer_df["mean_depth"] * infer_df["K_loess"]

    # Step 5: Calculate the average correction depth(mu_2) and depth per copy of diploid(a)
    diploid_df["corrected_depth"] = diploid_df["mean_depth"] * diploid_df["K_loess"]
    valid_diploid = diploid_df[diploid_df["mean_depth"] > 0]
    mu_2 = valid_diploid["corrected_depth"].mean()
    a = mu_2 / 2

    print(f"[INFO] mu_2 = {mu_2:.4f}, a = {a:.4f}")

    # Step 6: Copy number estimation
    infer_df["estimated_CN"] = infer_df["corrected_depth"] / a

    # Step 7: Output
    infer_df.to_csv(output_file, sep="\t", index=False)
    print(f"[DONE] Output {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GC content correction and copy number estimation.")
    parser.add_argument("-d", "--diploid", required=True, help="Diploid bin GC-depth file (e.g., *_bin_gc_depth.txt)")
    parser.add_argument("-i", "--input", required=True, help="Input window GC-depth file (e.g., *.chr*.gc.depth.bed)")
    parser.add_argument("-o", "--output", required=True, help="Output file name (e.g., *.gc.depth.cn.bed)")

    args = parser.parse_args()

    main(args.diploid, args.input, args.output)
