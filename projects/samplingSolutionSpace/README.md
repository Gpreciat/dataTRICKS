<p align="center">
  <img src="https://img.shields.io/badge/Metabolic%20Engineering-Expert-green.svg" alt="Metabolic Engineering">
  <img src="https://img.shields.io/badge/Software%20Development-Expert-green.svg" alt="Software Development">
  <img src="https://img.shields.io/badge/Research-Expert-green.svg" alt="Research">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Mathematical%20Optimization-Enthusiast-blueviolet.svg" alt="Mathematical Optimization">
</p>

# Sampling Repository

This repository provides a comprehensive collection of sampling algorithms for metabolic models, including methods such as CHRR (Coordinate Hit-and-Run with Rounding). These algorithms allow for uniform sampling of the flux space, enabling researchers to gain insights into metabolic network behavior.

## Usage

The sampling algorithms available in this repository can also be found in the [COBRA Toolbox](https://github.com/opencobra/cobratoolbox), a widely used platform for constraint-based modeling of metabolic networks.

### Identifying Significant Reactions

One of the key features of this repository is the `identifySignificantRxns` function, which leverages precomputed sampling data to identify the most crucial reactions within a genome-scale model. The process involves three key steps:

1. **Sampling the Flux Space:** Initiate the analysis by sampling the flux space.
2. **Calculating the Covariance Matrix:** Compute the covariance matrix based on the sampled data.
3. **Computing the Euclidean Norm:** Determine the Euclidean norm to quantify the significance of reactions.

This approach enables the identification of reactions for which constraining their bounds would have a substantial impact on the original flux space. The methodology was initially proposed by [Preciat et al.](https://www.biorxiv.org/content/10.1101/2021.06.30.450562v1) for determining which reactions to measure exometabolically to maximize the obtained information.

### Undersampling and Sampling

The repository provides tools for both undersampling and sampling of metabolic models. Undersampling with the CHRR algorithm can be performed using the `sampleCbModel` function. Sampling parameters should be carefully selected to ensure valid sampling distributions. Detailed instructions on parameter settings are provided in the documentation.

### Flux Variability Analysis

Additionally, this repository offers functionality for Flux Variability Analysis (FVA), which returns the minimum and maximum possible flux through every reaction in a model. It allows for a comprehensive exploration of the model's capabilities under various conditions.

## Getting Started

To get started with these sampling algorithms, make sure you have the required dependencies installed as specified in the [COBRA Toolbox](https://github.com/opencobra/cobratoolbox) installation guide. Initialize the COBRA Toolbox and verify that the pre-packaged LP and QP solvers are functional.

## Examples

Explore the provided examples to understand how to use these sampling algorithms effectively. From identifying significant reactions to performing FVA and sampling metabolic models, this repository offers a wide range of tools for your research needs.

