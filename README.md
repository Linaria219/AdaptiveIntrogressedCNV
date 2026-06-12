# AdaptiveIntrogressedCNV
[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![Linux](https://img.shields.io/badge/Platform-Linux-orange.svg)]()

# Adaptive Introgressed CNV Analysis Pipeline

## Overview
An end-to-end automated pipeline​ for large-scale genomic structural variation analysis.
The Challenge:​ Processing 2,000+ whole-genome samples​ to identify high-confidence Copy Number Variations (CNVs) while minimizing false positives inherent in single-algorithm approaches.
The Solution:​ A modular, multi-caller consensus framework written in Python/Bash, designed for HPC cluster environments. It automates data ingestion, parallel processing, statistical validation, and visualization.
Scale:​ Handles TB-scale datasets with robust error handling and logging.

## Pipeline Flowchart

```mermaid
graph TD
    classDef process fill:#e1f5fe,stroke:#01579b,stroke-width:2px;
    classDef data fill:#fff9c4,stroke:#fbc02d,stroke-width:1px;

    A([<b>Raw WGS/WES Data</b>]):::data
    
    subgraph S1 [Stage 1: Discovery]
        B[CNV Callers<br/><i>x4 Parallel Jobs on HPC</i>]
    end
    
    C{Consensus Filter<br/><b>>=2 Callers?</b>}:::process
    
    subgraph S2 [Stage 2: Quantification]
        D[dCGH Profiling<br/><i>Continuous CN Estimation</i>]
    end
    
    E[Stratification<br/><i>Classify by CN Patterns<br/>(Fixed/Polymorphic)</i>]:::process
    
    subgraph S3 [Stage 3: Evolutionary Insight]
        F[Positive Selection Scan<br/><i>iHS / Tajima's D</i>]
        G[Archaic Introgression<br/><i>100-SNV Sliding Window</i>]
    end
    
    H([<b>High-Confidence CNV Report</b>]):::data

    %% Links
    A --> B
    B --> C
    C -->|No (Discard)| X[Noise Filtered Out]
    C -->|Yes (High-Confidence BED)| D
    D --> E
    E --> F
    E --> G
    F --> H
    G --> H
```

## Tech Stack & Engineering Evidence
- **Languages:** Python (Pandas, NumPy, SciPy), R (ggplot2, data.table), Bash/Shell.
- **Environment:** Linux (CentOS/Ubuntu), HPC Clusters, SLURM Job Scheduler.
- **Data Scale:** 2000+ WGS samples(>=30x), ~200TB raw data processed.
- **Engineering Practices:**
  - Modular code structure for maintainability.
  - Batch automation scripts for cluster computing.
  - Standardized output formats (BED, VCF, CSV).

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
- Auxiliary scripts used for data parsing, format conversion, and information collection.

### 7. Visualization (`R/`)
- Helper scripts for visualization.

## Tech Stack & Evidence
- **Languages**: Python, R, Bash (Linux environment).
- **Workflow Management**: Modular shell scripting for batch processing.

