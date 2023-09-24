<p align="center">
  <img src="https://img.shields.io/badge/Metabolic%20Engineering-Expert-green.svg" alt="Metabolic Engineering">
  <img src="https://img.shields.io/badge/Software%20Development-Expert-green.svg" alt="Software Development">
  <img src="https://img.shields.io/badge/Research-Expert-green.svg" alt="Research">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Mathematical%20Optimization-Enthusiast-blueviolet.svg" alt="Mathematical Optimization">
</p>

# Dopaminergic Neuron GEM (iDopaNeuro): Omics-Powered Repository

<p align="center">
  <img src="data/neuron.png" alt="Neuron">
</p>

## Introduction

This repository contains the code and resources related to the research project aimed at gaining a deeper understanding of the metabolic mechanisms underlying Parkinson's Disease. The work is currently under review in *Nature Communications in Biology* (Link to [BioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.30.450562v2)).

## Project Overview

The objective is to design a Genome-Scale Metabolic Model (GEM) specific to dopaminergic neurons. What sets our methodology apart is the integration of omics data, including transcriptomics, metabolomics, and bibliomics, into the GEM construction process. It was also developed a custom function called [XomicsToModel](https://github.com/Gpreciat/dataTRICKS/tree/main/projects/omicsDataIntegration) that plays a crucial role in generating omics-integrated GEMs.

## Methodology

The methodology involves the following key steps:

1. **In Vitro Differentiation**: Human neuroepithelial stem cells were differentiated into midbrain-specific dopaminergic neuronal cultures in vitro.

2. **Data Generation**: Transcriptomic and targeted exometabolomic data were generated from fresh and spent media samples.

3. **Literature Curation**: Cell-type-specific data were derived from manual curation of the literature on dopaminergic neuronal metabolism.

4. **Model Integration**: Omics data were integrated with a generic metabolic model derived from a comprehensive reconstruction of human metabolism.

5. **Model Generation**: An ensemble of candidate dopaminergic neuronal metabolic models were generated, considering technical parameters and mathematical modeling approaches.

6. **Model Selection**: Models with the highest predictive fidelity were identified, giving preference to either in vitro experimental data (iDopaNeuro condition specific; iDopaNeuroC) or dopaminergic neuronal literature curation (iDopaNeuro cell-type specific; iDopaNeuroCT).

7. **Validation**: The predictive fidelity of the iDopaNeuro models was validated by comparison with independent exometabolomic data generated on perturbations to normal dopaminergic neuronal metabolism in vitro.

8. **Prospective Experiments**: Finally, the iDopaNeuro models were used prospectively to design exometabolomic experiments to generate new constraints on the variables currently [most uncertain](https://github.com/Gpreciat/dataTRICKS/tree/main/projects/samplingSolutionSpace) in the model.

## Conclusion

This multidisciplinary project aims to shed light on the metabolic mechanisms of Parkinson's Disease using cutting-edge computational modelling and data integration techniques. This work is a testament to the power of combining diverse sources of data to gain insights into complex biological systems.

Please explore the code and resources in this repository to learn more about the research.

*Note: The preprint of this work is currently under review in Nature Communications in Biology, and the link provided is to the preprint version.*
