# Protein-Ligand Proximity Analysis

This Python script performs a proximity analysis between specified residues in a protein and the center of a ligand. It determines if the residues are within various cutoff radii from the ligand's geometric center. The analysis results are outputted to a CSV file, which includes whether each residue is within the specified cutoff radii.

## Features

- Parses PDB files for protein and ligand structures.
- Calculates the geometric center of the ligand.
- Checks proximity of specified residues to the ligand center for various cutoff radii.
- Outputs the results to a CSV file with detailed information about each residue's proximity to the ligand.

## Prerequisites

Before running this script, ensure you have Python installed on your system along with the following packages:
- `numpy`
- `biopython`
- `pyyaml`

You can install these packages using pip:

```bash
pip install numpy biopython pyyaml
