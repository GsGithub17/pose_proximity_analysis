#Pose Proximity Analysis
'''A Python function that reads in PDB structures for a protein and a ligand, calculates the distance between the center of the ligand and specified residues in the protein, and then determines whether this distance is within various cutoff radii and returns Boolean values. It’s mainly based on Biopython to use the readily available PDBParser function for both the target and ligand pdb structures.'''
#Author Gemechis Degaga (Ph.D, email:gddegaga@gmail.com)
import numpy as np
from Bio.PDB import PDBParser, Selection
import csv
import yaml
import os
import re  # For parsing the custom residue ID format

#function to calculate center of a given ligand or residues.
def calculate_center(atoms):
    coords = [atom.get_coord() for atom in atoms]
    return np.mean(coords, axis=0)

#function toParse the custom key_residues formatted like 'A_Met111'.
def parse_residue_id(residue_id):
    pattern = re.compile(r'([A-Za-z])_([A-Za-z]{3})(\d+)')
    match = pattern.match(residue_id)
    #the following conditional block checks whiter a list of proper key_residue list are provided.
    if match:
        chain = match.group(1)
        resname = match.group(2)
        resnum = int(match.group(3))
        return chain, resname, resnum
    else:
        raise ValueError("Residue ID format is incorrect")

#funcition for actual proximity analysis
def proximity_analysis(config_file):
    #getting current working directory, input path paths for the target and ligand are relative to this working directory.
    current_directory = os.getcwd()
    #reading in the inputs for the config.yaml file
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    
    protein_pdb_path = os.path.join(current_directory, config['protein_pdb_path'])
    ligand_pdb_path = os.path.join(current_directory, config['ligand_pdb_path'])
    key_residues = config['key_residues']
    cutoff_radii = config['cutoff_radii']

    parser = PDBParser(QUIET=True)
    structure_protein = parser.get_structure('Protein', protein_pdb_path)
    structure_ligand = parser.get_structure('Ligand', ligand_pdb_path)
    
    ligand_atoms = Selection.unfold_entities(structure_ligand, 'A')  # All atoms in the ligand
    ligand_center = calculate_center(ligand_atoms)
    
    model_protein = structure_protein[0]  # Assuming the protein has only one model

    results = []

    for key_residue in key_residues:
        chain_id, resname, resnum = parse_residue_id(key_residue)
        chain = model_protein[chain_id]
        residue = next((res for res in chain if res.get_id()[1] == resnum and res.get_resname() == resname), None)
        if residue:
            residue_atoms = residue.get_unpacked_list()
            residue_center = calculate_center(residue_atoms)
            distances = [np.linalg.norm(residue_center - ligand_center) <= r for r in cutoff_radii]
            results.append([ligand_pdb_path.split('/')[-1], chain_id, residue.get_resname(), residue.get_id()[1]] + distances)
        else:
            print(f"Residue {resname}{resnum} not found in chain {chain_id}")
    
    # Write results to CSV
    header = ['Ligand', 'Chain', 'Residue_Name', 'Residue_ID'] + [f'within_{r}A' for r in cutoff_radii]
    with open('proximity_analysis_results.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)
        writer.writerows(results)

# Example usage
proximity_analysis('config.yaml')

