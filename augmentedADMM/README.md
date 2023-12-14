# Standard ADMM vs. Augmented ADMM Comparison

This R package provides a comprehensive comparison between the Standard Alternating Direction Method of Multipliers (ADMM) and the Augmented ADMM. It includes a series of scripts for testing and evaluating the performance of these algorithms in various scenarios, including standard settings and graph-structured data.

## Overview

The package contains scripts that compare Standard ADMM and Augmented ADMM in terms of runtime and precision. Each script generates a data set, applies both ADMM algorithms, and plots the results for easy visualization and analysis. The comparison covers both standard data matrices and graph-structured data, providing a broad view of the performance characteristics of these algorithms.

## Guide to Using the Package

A detailed guide is provided in the form of an RMarkdown file (`guide.Rmd`). This file includes instructions for installing the package.

## Getting Started

To run these scripts, clone the repository and ensure you have R and the necessary libraries installed. Each script can be executed independently to generate the corresponding plots.

## Test Files

### Runtime Comparison Scripts

These scripts compare the runtime efficiency of Standard ADMM and Augmented ADMM. 

- **`std_vs_aug_admm_runtime_comparison.R`**: 
  - **Purpose**: Compares the runtime performance on standard data sets.
  - **Output**: `std_vs_aug_admm_runtime_comparison.pdf`, a plot visualizing the runtime comparison.

- **`std_vs_aug_for_graph_runtime_comparison.R`**: 
  - **Purpose**: Focuses on runtime comparison using graph-structured data, based on the graph generating procedure from Zhu, Y. (2017).
  - **Output**: `std_vs_aug_admm_for_graph_runtime_comparison.pdf`, a plot illustrating runtime comparison for graph data.

### Precision Comparison Scripts

These scripts evaluate the precision of the two ADMM methods.

- **`std_vs_aug_admm_precision_comparison.R`**: 
  - **Purpose**: Assesses the precision of Standard and Augmented ADMM on standard data sets.
  - **Output**: `std_vs_aug_admm_precision_comparison.pdf`, a plot showing the precision comparison.

- **`std_vs_aug_for_graph_precision_comparison.R`**: 
  - **Purpose**: Analyzes the precision of the algorithms specifically for graph-structured data, using data generation techniques from Zhu, Y. (2017).
  - **Output**: `std_vs_aug_admm_for_graph_precision_comparison.pdf`, a plot detailing the precision comparison for graph data.



## References

- Zhu, Y. (2017). An Augmented ADMM Algorithm With Application to the Generalized Lasso Problem. Journal of Computational and Graphical Statistics.
