<p align="center">
  <img src="https://img.shields.io/badge/Metabolic%20Engineering-Expert-green.svg" alt="Metabolic Engineering">
  <img src="https://img.shields.io/badge/Software%20Development-Expert-green.svg" alt="Software Development">
  <img src="https://img.shields.io/badge/Research-Expert-green.svg" alt="Research">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Mathematical%20Optimization-Enthusiast-blueviolet.svg" alt="Mathematical Optimization">
</p>

# XomicsToModel - Omics data integration in genome-scale models

## Introduction

The `XomicsToModel` function is a powerful tool within the COBRA Toolbox designed to bridge the gap between multi-omics data and genome-scale metabolic models. This function serves as a transformative bridge, allowing researchers to translate various omics datasets, such as transcriptomics, proteomics, metabolomics, and literature-based "bibliomics," into context-specific metabolic models for particular cell types or conditions. This methodology was utilised in the study by *Preciat et al.* to construct a genome-scale model of cell dopaminergic neurons derived from pluripotent stem cells ([iDopaNeuro](https://github.com/Gpreciat/dataTRICKS/tree/main/projects/dopaminergicNetworkGem)). The omics information employed encompassed bibliomic, metabolomic, and transcriptomic data.

## Key Features

- **Multi-Omics Integration:** `XomicsToModel` seamlessly integrates diverse omics data types, including transcriptomics, proteomics, metabolomics, and bibliomics, to derive context-specific metabolic models.

- **Contextualization:** It generates metabolic models tailored to specific cell types or conditions, providing insights into the metabolic behaviors that are relevant to the given context.

- **Flexible Input:** Accepts a generic COBRA model and context-specific omics data as input, making it adaptable to a wide range of research scenarios.

- **Parameter Customization:** Offers a variety of parameters for fine-tuning the model generation process, allowing researchers to control and optimize their analyses.

## How to Use

To utilize the `XomicsToModel` function, follow these simple steps:

```matlab
[model, modelGenerationReport] = XomicsToModel(genericModel, specificData, param);
```

Where:
- `genericModel`: A generic COBRA model.
- `specificData`: A structure containing context-specific omics data.
- `param`: A structure containing various parameters for customization.

## Getting Started

You can also find scripts and functions related to this project in the [COBRA Toolbox](https://github.com/opencobra/cobratoolbox). For a detailed guide on using the `XomicsToModel` function, please refer to the [COBRA Toolbox Documentation](https://github.com/opencobra/cobratoolbox/blob/master/docs/source/notes/COBRAModelFields.md).

## Publication (Manuscript Under Review)

For a comprehensive understanding of the `XomicsToModel` function and its applications, please refer to our publication (manuscript under review). We are currently in the review process, and the manuscript will be updated to the corresponding journal once it is accepted.

[Link to Manuscript (Under Review)](https://www.biorxiv.org/content/10.1101/2021.11.08.467803v1)

## License

This project is dual-licensed under two licenses: [dataTRICKS](https://github.com/Gpreciat/dataTRICKS/blob/main/LICENSE.txt) (since I authored the code) and the [COBRA Toolbox](https://github.com/opencobra/cobratoolbox/tree/master/) license (as the project was further developed there). Users are encouraged to review and comply with the terms of both licenses.

For details on the licenses, please refer to the respective repositories.

## Contribution and Contact

We welcome contributions and feedback from the research community. If you have any questions, suggestions, or would like to collaborate, please feel free to contact us through the information provided below:

**German Preciat, PhD**
- LinkedIn: [@gpreciat](https://www.linkedin.com/in/gpreciat/)
- Email: gapreciat@gmail.com

Thank you for choosing `XomicsToModel`. Together, let's unravel the mysteries of metabolism! 🧬🔬📊
