
# Chemoinformatics Database Generation

## Introduction

The Chemoinformatics Database Generation project offers a powerful tool for generating a chemoinformatics database from a genomic-scale model. By leveraging metabolite identifiers found in GEMs (Genome-scale Metabolic Models), this project enables the retrieval of molecular structures for metabolites from various databases or in different formats (e.g., MOL, RXN, SMILES, InChI, RInChI and InChIKey). 

## Features

- **Metabolite Identification:** Obtain molecular structures for metabolites.
- **Structural Comparison:** Compare molecular structures based on their InChI representation.
- **In-Depth Comparison:** Includes chemical formula, charge, stereochemistry, and similarity with other databases.
- **Standardization:** Check if the correct pH.

## Usage

1. **Metabolite Retrieval:** Using metabolite identifiers from GEMs, retrieve molecular structures from databases.

2. **Structural Comparison:** Compare retrieved structures to find the most suitable one for your genomic-scale model.

3. **Reaction Generation:** Generate atomically balanced reactions along with files in RXN, RInChI, and SMILES formats.

## Integration with COBRA Toolbox

You can also find scripts and functions related to this project in the COBRA Toolbox. For more detailed results using different genomic-scale models, please refer to [the GitHub repository here](https://github.com/opencobra/ctf).

## Applications

- **Atom Transition Networks at Genome Scale:** Understand atomic transitions in a genomic context.
- **Identification of Conserved Moieties:** Discover conserved biochemical substructures.
- **Fluxomics Experiment Design:** Design experiments to analyze cellular fluxes.
- **Functional Objective Design:** Develop functional objectives for metabolic models.
- **Identification of Enthalpy and Bond Changes:** Investigate thermodynamic properties at a genomic scale.

## Getting Started

All project code and functions are available in the [COBRA Toolbox](https://github.com/opencobra/cobratoolbox/tree/master/). To get started, clone the COBRA Toolbox repository and explore the `chemoInformatics` folder.

## License

This project is dual-licensed under two licenses: [dataTRICKS](https://github.com/Gpreciat/dataTRICKS/blob/main/LICENSE.txt) (since I authored the code) and the [COBRA Toolbox](https://github.com/opencobra/cobratoolbox/tree/master/) license (as the project was further developed there). Users are encouraged to review and comply with the terms of both licenses.

For details on the licenses, please refer to the respective repositories.

Feel free to contribute to this project and leverage its capabilities to enhance your genomic-scale metabolic modeling tasks.

Happy coding!