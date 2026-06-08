# AdaptiveIntrogressedCNV
# Adaptive Introgressed CNV Analysis Pipeline

## Overview
This repository contains an end-to-end bioinformatics pipeline for identifying **high-confidence Copy Number Variations (CNVs)** and investigating their evolutionary significance, specifically focusing on **positive selection** and **archaic introgression**.

## Pipeline Architecture
The workflow is modularized into six stages:

### 1. CNV Discovery & Consensus (`CNVcallers/`)
- **Input**: Raw WGS/WES data.
- **Method**: Runs 4 independent CNV calling algorithms.
- **Filtering**: Implements an iterative merging strategy. Only CNVs detected by **≥2 callers** are retained as the high-confidence dataset for downstream analysis.
- **Output**: Merged BED files of reliable CNV loci.

### 2. Quantitative CNV Profiling (`dCGH_pipeline/`)
- **Method**: Applies density-based Copy Number Variation analysis using Genomic Hybridization (dCGH).
- **Goal**: Estimates the actual copy number (CN) value for each CNV region across all samples, moving beyond binary "gain/loss" calls.

### 3. CNV Stratification (`Stratification/`)
- **Analysis**: Classifies CNVs based on their CN distribution patterns (e.g., fixed, polymorphic, or population-specific).

### 4. Positive Selection Scan (`Adaptation/`)
- **Method**: Identifies Single Nucleotide Variants (SNVs) under positive natural selection using population genetics metrics.
- **Integration**: Tests for overlap between selected SNVs and CNV regions to infer functional constraint.

### 5. Archaic Introgression Analysis (`Introgression/`)
- **Method**: Scans 100-SNV sliding windows for signals of Neanderthal/Denisovan ancestry.
- **Goal**: Determines if high-confidence CNVs reside within regions of archaic introgression.

### 6. Utilities (`Tools/`)
- Helper scripts for data parsing, format conversion, and visualization.

## Tech Stack & Evidence
- **Languages**: Python, R, Bash (Linux environment).
- **Workflow Management**: Modular shell scripting for batch processing.

