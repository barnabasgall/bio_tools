import glob
import pandas as pd
import networkx as nx
import itertools
from Bio.PDB import Superimposer, PDBParser
from Bio.PDB import PDBIO, Superimposer, PDBParser
import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDWriter
from MDAnalysis import Universe
import tempfile
import argparse
from tqdm import tqdm


def greedy_hamiltonian_path(G):
    current_node = list(G.nodes)[0]  # start at a random node
    path = [current_node]
    unvisited = set(G.nodes)
    unvisited.remove(current_node)

    while unvisited:
        next_node = min(unvisited, key=lambda x: G[current_node][x]['rmsd'])
        unvisited.remove(next_node)
        path.append(next_node)
        current_node = next_node

    return path

def calculate_c_alpha_rmsd(structure_1, structure_2):
    super_imposer = Superimposer()

    atoms_1 = [atom for atom in structure_1.get_atoms() if atom.get_name() == 'CA']
    atoms_2 = [atom for atom in structure_2.get_atoms() if atom.get_name() == 'CA']

    if len(atoms_1) != len(atoms_2):
        raise ValueError("Structures do not have the same number of C-alpha atoms!")

    super_imposer.set_atoms(atoms_1, atoms_2)

    return super_imposer.rms

def calculate_c_alpha_rmsd_from_files(file_1, file_2):
    parser = PDBParser()
    structure_1 = parser.get_structure('structure_1', file_1)
    structure_2 = parser.get_structure('structure_2', file_2)

    return calculate_c_alpha_rmsd(structure_1, structure_2)

def get_atoms(filename, get_struc = False):
    parser = PDBParser()
    structure = parser.get_structure("temp", filename)
    atoms = list(structure.get_atoms())
    if get_struc:
        return atoms, structure
    else:
        return atoms

def calc_avg_bfactor(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)

    b_factors = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    b_factors.append(atom.get_bfactor())
    
    avg_bfactor = sum(b_factors) / len(b_factors)
    
    return avg_bfactor

def make_dcd(folder, output_name, filter = True, cutoff=85, file_expression = '*_rel*.pdb'):
    files = glob.glob(folder + '/' + file_expression)

    if filter:
        files = [file for file in files if calc_avg_bfactor(file) > cutoff]

    df = pd.DataFrame(columns=['file_1', 'file_2', 'rmsd'])

    # make RMSD

    for file_1, file_2 in tqdm(itertools.combinations(files, 2), 'calculating pairwise RMSD'):
        rmsd = calculate_c_alpha_rmsd_from_files(file_1, file_2)
        df = pd.concat([df,pd.DataFrame({'file_1': [file_1], 'file_2': [file_2], 'rmsd': [float(rmsd)]})])
    
    df.to_csv(output_name+'.csv')


    G = nx.from_pandas_edgelist(df, 'file_1', 'file_2', 'rmsd')
    path = greedy_hamiltonian_path(G)

    io = PDBIO()

    # Create an empty Universe
    u = mda.Universe.empty(0)

    # Initialize Superimposer
    super_imposer = Superimposer()

    # Get atoms of the first structure
    ref_atoms, structure = get_atoms(path[0], get_struc = True)
    # save pdb file
    io.set_structure(structure) 
    io.save(output_name+'.pdb')

    # Initialize DCD writer
    with DCDWriter(output_name+'.dcd', n_atoms=len(ref_atoms)) as W:
        for pdb_file in tqdm(path,'writing dcd'):
            # check b-factor
            avg_bfactor = calc_avg_bfactor(pdb_file)
            if avg_bfactor < 85:
                continue
            else:
                # Parse the structure
                structure = PDBParser().get_structure('temp', pdb_file)
                atoms = list(structure.get_atoms())

                # Superimpose the structure on the reference structure
                if len(ref_atoms) != len(atoms):
                    print("Reference and target atom lists have different sizes. Skipping superimposition for this file.")
                    continue
                super_imposer.set_atoms(ref_atoms, atoms)
                super_imposer.apply(atoms)

                # Create an MDAnalysis Universe with the superimposed structure
                u = Universe.empty(n_atoms=len(atoms), trajectory=True)
                u.atoms.positions = [atom.get_coord() for atom in atoms]

                # Write the universe to the DCD
                W.write(u.atoms)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', type=str, help='folder containing pdbs')
    parser.add_argument('--output_name', type=str, help='output name')
    parser.add_argument('--filter', type=bool, help='filter by b-factor', default=True)
    parser.add_argument('--cutoff', type=int, help='b-factor cutoff', default=85)
    parser.add_argument('--file_expression', type=str, help='file expression', default='*_rel*.pdb')
    args = parser.parse_args()
    make_dcd(args.folder, args.output_name)