from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB import PDBIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.PDB.Polypeptide import three_to_one
from Bio.Alphabet import generic_protein
from itertools import permutations

import numpy as np
import math

from overlapy.needleman import Needleman


class Superimposer:
    """
    This class is used to perform superimposition between 2 PDBs.

    It follows the extract method to method object design pattern. You enter all the input parameters,
    then you call superimpose() and you collect the results. You will have to create a superimposition class
    for each superimposition. YOU CAN'T REUSE INSTANCES.

    Sample usage:

        >>> si = Superimposer()
        >>> si.parse("pdb_1.pdb", "pdb_2.pdb")
        >>> si.select_chains(["A", "B"], ["A", "B"]) # optional
        >>> si.superimpose() # superimpose method
        >>> rmsd = si.rmsd
        >>> rotation_matrix = si.rotation_matrix
        >>> ms_alignment = si.get_multiple_sequence_alignment()
        >>> si.save_superimposed_pdb("superimposed.pdb")
    """

    def __init__(self):
        self.rotation_matrix = None
        self.rmsd = None
        self.alignment = None
        self.already_executed = False

    # INPUT PARAMETERS

    def parse(self, *pdb_filenames):
        """
        REQUIRED. Adds the protein PDB files. You can specify as many as you want, but only two will be used for the superimposition.
        """
        self.proteins = [] # reset proteins to an empty array
        parser = PDBParser(QUIET=True)
        for filename in pdb_filenames:
            # use file name as PDB id
            pdb_id = self.__get_pdb_id_from_filename(filename)
            # get PDB contents
            self.proteins.append(parser.get_structure(pdb_id, filename))

    def select_chains(self, *chain_lists):
        """
        OPTIONAL. Restrains the superimposition to the given chains in each protein. Chain lists have to be specified as lists.

        Example:
            >>> superimposer.select_chains(['A'], ['A', 'B'])
        """
        self.used_chains = chain_lists
        for i,chain_list in enumerate(chain_lists):
            # check that the list is not empty
            if chain_list:
                self.__remove_unused_chains(i, chain_list)

    # ACTION

    def superimpose(self):
        if not self.already_executed:
            # center protein atoms
            for id,protein in enumerate(self.proteins):
                self.__move_atoms_to_center(id)
            # perform needleman
            self.alignment = self.__needleman()
            # get rotation matrix
            self.rotation_matrix = self.__compute_rotation_matrix(self.alignment)

            if self.__valid_alignment():
                # rotate second protein
                self.__rotate_protein(1)
                # calculate RMSD
                self.rmsd = self.__compute_rmsd()
                # only run once superimpose
            self.already_executed = True

    # OUTPUT PARAMETERS

    def save_superimposed_pdb(self, out_filename):
        """
        Saves the superimposed PDB in the given output filename.
        """
        if self.__valid_alignment():
            superimposed_pdb = self.__create_superimposed_pdb()

            # save it to a file
            io = PDBIO()
            io.set_structure(superimposed_pdb)
            io.save(out_filename)

    def get_multiple_sequence_alignment(self):
        """Returns a MultipleSeqAlignment object with the alignment given by Needleman & Wunsch."""
        def create_seqrecord(sequence, name):
            sequence_string = ""
            for aa in sequence:
                if aa is None:
                    symbol = "-"
                else:
                    try:
                        symbol = three_to_one(aa.get_resname())
                    except:
                        symbol = "?"
                sequence_string += symbol
            return SeqRecord(Seq(sequence_string, generic_protein), id=name)

        first_sequence = create_seqrecord(self.alignment[0], "First")
        second_sequence = create_seqrecord(self.alignment[1], "Second")
        alignment = MultipleSeqAlignment([first_sequence, second_sequence])
        return alignment

    # PRIVATE

    def __get_pdb_id_from_filename(self, pdb_filename):
        return pdb_filename.split('/')[-1].split('.')[0]

    def __move_atoms_to_center(self, protein_id):
        atoms = list(self.proteins[protein_id].get_atoms())
        centroid = self.__compute_centroid(atoms)
        # centroid = self.__compute_centroid((atom for atom in atoms if atom.id == 'CA'))
        for atom in atoms:
            atom.set_coord(atom.get_coord() - centroid)

    def __compute_centroid(self, atoms):
        """
        Given a list of atoms, return the geometrical center of all of them.
        This is, the average of each of the coordinates.
        """
        matrix = self.__atom_coords_to_matrix(atoms)
        centroid = []
        for column in range(matrix.shape[1]):
            centroid.append(np.mean(matrix.T[column]))
        return np.array(centroid)

    def __atom_coords_to_matrix(self, atoms):
        return np.matrix([atom.get_coord() for atom in atoms])

    def __compute_score_matrix(self, first_sequence, second_sequence):
        """
        Given two sequences, return a matrix of scores. If the
        structures have n and m alpha carbons respectively, the matrix
        of scores will be of dimensions n by m.
        The cell matrix[i][j] will contain the score of the ith alpha carbon
        from structure 1 and jth alpha carbon from structure2.

        """
        alphas1 = self.__get_alpha_carbons(first_sequence)
        alphas2 = self.__get_alpha_carbons(second_sequence)

        score_matrix = np.empty((len(alphas1)-10+1, len(alphas2)-10+1))

        # iterate through all alpha carbons, but start at 5 and end at -4
        # (we need an environment around each of them)
        for i, alpha1 in enumerate(alphas1[5:-4]):
            distances_1 = self.__get_environment_distances(alphas1, i, 5)
            for j, alpha2 in enumerate(alphas2[5:-4]):
                distances_2 = self.__get_environment_distances(alphas2, j, 5)
                score_matrix[i][j] = np.abs(distances_1 - distances_2).sum()

        return score_matrix

    def __get_environment_distances(self, atoms, index, radius):
        """
        Given a list of atoms, an index and a radius (in number of atoms),
        return a numpy array of distances from the atom at the given index to
        each of the atoms between [index - radius, index + radius + 1].

        The distance of the atom at 'index' to itself is not included.

        """
        distances = np.empty((radius*2))
        center_atom = atoms[index]
        # local list removing the atom at position index
        local_atoms = []
        for i in range(index-radius, index+radius+1):
            if i != index:
                local_atoms.append(atoms[i])
        for i, atom in enumerate(local_atoms):
            distances[i] = atom - center_atom
        return distances

    def __get_aminoacids(self, protein_id, chains):
        residues = []
        protein = self.proteins[protein_id][0]
        for chain in chains:
            residues.extend([residue for residue in protein[chain].get_residues() if 'CA' in residue])
        return residues

    def __get_alpha_carbons(self, aa_sequence):
        return [residue['CA'] for residue in aa_sequence]

    def __needleman(self):
        # first and last 5 aminoacids are not considered for the alignment
        chains_first = self.__get_all_chains(0)
        chains_second = self.__get_all_chains(1)

        # all permutations of the chains
        perms_first = list(permutations(chains_first))
        perms_second = list(permutations(chains_second))

        best_score = None
        # aling all the permutations of the chains and return the best alignment
        for perm_first in perms_first:
            for perm_second in perms_second:
                # ordered residues, according to the current permutations
                aminoacids1 = self.__get_aminoacids(0, perm_first)[4:-5]
                aminoacids2 = self.__get_aminoacids(1, perm_second)[4:-5]

                score_matrix = self.__compute_score_matrix(aminoacids1, aminoacids2)

                needleman = Needleman(aminoacids1, aminoacids2, score_matrix)
                needleman.align()
                score = needleman.score
                if best_score is None or score < best_score:
                    best_score = score
                    best_alignment = needleman.get_alignment()

        return best_alignment


    def __compute_rotation_matrix(self, alignment):
        """
        Given the structural alignment, compute the rotation matrix which
        minimizes the RMSD among alpha carbons.
        The alignment must be passed as a tuple of two list, each corresponding
        to one of the sequences. Each list is composed of aligned residues
        (objects of type residue) or None to represent a gap.

        """
        aminoacids1 = []
        aminoacids2 = []
        for i in range(len(alignment[0])):
            if alignment[0][i] and alignment[1][i]:
                aminoacids1.append(alignment[0][i])
                aminoacids2.append(alignment[1][i])

        # represent each list of aminoacids as a matrix of coordinates of the
        # alpha carbons.
        matrix1 = self.__atom_coords_to_matrix((aa['CA'] for aa in aminoacids1))
        matrix2 = self.__atom_coords_to_matrix((aa['CA'] for aa in aminoacids2))

        if matrix1.any() and matrix2.any():
            # build matrix of covariances: transpose of the first, by the second
            matrix_covariances = matrix1.T * matrix2

            # SVD decomposition of the covariance matrix. Wt = transpose of W
            V, S, Wt = np.linalg.svd(matrix_covariances)

            det_sign = 1 if np.linalg.det(Wt.T * V.T) > 0 else -1

            # rotation matrix = W * transpose of T
            return Wt.T * np.array([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, det_sign]]) * V.T
        else:
            return None

    def __create_superimposed_pdb(self):
        def fill_in_chain(chain, protein_id, rotation_matrix = None):
            for index,residue in enumerate(self.proteins[protein_id].get_residues()):
                residue.id = (residue.id[0], index, residue.id[2])
                chain.add(residue)

        merged_model = Model(0)
        chain_a = Chain('A')
        chain_b = Chain('B')

        fill_in_chain(chain_a, 0)
        fill_in_chain(chain_b, 1)

        merged_model.add(chain_a)
        merged_model.add(chain_b)

        return merged_model

    def __rotate_protein(self, protein_id):
        for atom in self.proteins[protein_id].get_atoms():
            coordinates = atom.get_coord()
            coordinates = np.dot(coordinates, self.rotation_matrix)
            # convert it again to a vector
            coordinates = np.squeeze(np.asarray(coordinates))
            atom.set_coord(coordinates)

    def __compute_rmsd(self):
        square_dists = [(res1['CA'] - res2['CA'])**2 for res1, res2 in zip(*self.alignment) if res1 and res2]
        return math.sqrt(sum(square_dists) / len(square_dists))

    def __remove_unused_chains(self, protein_id, used_chains):
        # we have to iterate through all models and ask if the chain exists
        for model in self.proteins[protein_id]:
            for model_chain in model:
                if model_chain.id not in used_chains:
                    model.detach_child(model_chain.id)

    def __get_all_chains(self, protein_id):
        return [chain.id for chain in self.proteins[protein_id].get_chains()]

    def __valid_alignment(self):
        return self.rotation_matrix is not None and self.rotation_matrix.any()
