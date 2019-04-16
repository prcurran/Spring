"""
The following script calculates:
    - natural product likeness score
    - synthetic accessibility score
    - plots scores as 2D KDE
    - commandline usage

This has been created using RDKit contributed scripts

P R CURRAN
"""
# calculation of synthetic accessibility score as described in:
#
# Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions
# Peter Ertl and Ansgar Schuffenhauer
# Journal of Cheminformatics 1:8 (2009)
# http://www.jcheminf.com/content/1/1/8
#
# several small modifications to the original paper are included
# particularly slightly different formula for marocyclic penalty
# and taking into account also molecule symmetry (fingerprint density)
#
# for a set of 10k diverse molecules the agreement between the original method
# as implemented in PipelinePilot and this implementation is r2 = 0.97
#
# peter ertl & greg landrum, september 2013
#
#
# calculation of natural product-likeness as described in:
#
# Natural Product-likeness Score and Its Application for Prioritization of Compound Libraries
# Peter Ertl, Silvio Roggo, and Ansgar Schuffenhauer
# Journal of Chemical Information and Modeling, 48, 68-74 (2008)
# http://pubs.acs.org/doi/abs/10.1021/ci700286x
#
# for the training of this model only openly available data have been used
# ~50,000 natural products collected from various open databases
# ~1 million drug-like molecules from ZINC as a "non-NP background"
#
# peter ertl, august 2015
#
from __future__ import division, print_function

import argparse
import gzip
import math
import os
import pickle
from collections import namedtuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

import seaborn as sns
import matplotlib.pyplot as plt


def fingerprint(mol):
    fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
    return fp.GetNonzeroElements()


class NaturalProductScorer(object):
    """
    calculates the natural product score for a given molecule
    """
    def __init__(self, npmodel):
        self.npmodel = npmodel

    def score(self, mol):
        """Next to the NP Likeness Score, this function outputs a confidence value
        between 0..1 that descibes how many fragments of the tested molecule
        were found in the model data set (1: all fragments were found).

        Returns namedtuple NPLikeness(nplikeness, confidence)"""

        if mol is None:
            raise ValueError('invalid molecule')
        bits = fingerprint(mol)

        # calculating the score
        score = 0.0
        bits_found = 0
        for bit in bits:
            if bit in self.npmodel:
                bits_found += 1
                score += self.npmodel[bit]

        score /= float(mol.GetNumAtoms())
        confidence = float(bits_found / len(bits))

        # preventing score explosion for exotic molecules
        if score > 4:
            score = 4. + math.log10(score - 4. + 1.)
        elif score < -4:
            score = -4. - math.log10(-4. - score + 1.)
        np_likeness = namedtuple("NPLikeness", "nplikeness,confidence")
        #return np_likeness(score, confidence)
        return np_likeness(score, confidence).nplikeness


class SyntheticAccesibility(object):
    """
    calculates the synthetic accessibility score for a given molecule
    """
    def __init__(self, fpscores):
        self.fragment_scores = fpscores

    def score(self, mol):
        """
        given an RDKit molecule a score is returned

        :param mol:
        :return:
        """

        # fragment score
        fps = fingerprint(mol)
        freqs, scores = zip(*[[v, self.fragment_scores.get(bitId, -4) * v] for bitId, v in fps.items()])
        fragment_score = sum(scores)/sum(freqs)

        # features score
        n_atoms = mol.GetNumAtoms()
        n_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        ri = mol.GetRingInfo()
        n_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        n_bridgeheads = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        n_macrocycles = len([x for x in ri.AtomRings() if len(x) > 8])


        size_penalty = n_atoms ** 1.005 - n_atoms
        stereo_penalty = math.log10(n_chiral_centers + 1)
        spiro_penalty = math.log10(n_spiro + 1)
        bridge_penalty = math.log10(n_bridgeheads + 1)
        macrocycle_penalty = 0.

        # This differs from the paper, which defines:
        #       macrocycle_penalty = math.log10(n_macrocycles+1)
        #
        # This form generates better results when 2 or more macrocycles are present

        if n_macrocycles > 0:
            macrocycle_penalty = math.log10(2)

        feature_score = 0. - size_penalty - stereo_penalty - spiro_penalty - bridge_penalty - macrocycle_penalty

        # correction for the fingerprint density
        # not in the original publication, added in version 1.1
        # to make highly symmetrical molecules easier to synthetise

        correction = 0.
        if n_atoms > len(fps):
            correction = math.log(float(n_atoms) / len(fps)) * .5

        sascore = fragment_score + feature_score + correction

        # need to transform "raw" value into scale between 1 and 10
        min = -4.0
        max = 2.5
        sascore = 11. - (sascore - min + 1) / (max - min) * 9.
        # smooth the 10-end
        if sascore > 8.:
            sascore = 8. + math.log(sascore + 1. - 9.)
        if sascore > 10.:
            sascore = 10.0
        elif sascore < 1.:
            sascore = 1.0

        return sascore


class Organiser(argparse.ArgumentParser):
    """
    class organising the natural product likeness and synthetic accessibility calculations
    """
    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)

        self.add_argument(
            'inputs',
            help='Path to input directory. The directory should contain all molecule sets and no other data!'
        )

        self.add_argument(
            'out_dir',
            help='path to output figure'
        )

        self.natural_product_scorer = NaturalProductScorer(self.read_np_model())
        self.synthetic_accessibility = SyntheticAccesibility(self.read_fragment_scores())

        self.args = self.parse_args()

        self.supported_file_formats = [".sdf"]
        self.input_files = [x for x in os.listdir(self.args.inputs)
                            if os.path.splitext(x)[1] in self.supported_file_formats]

        self.molecule_sets = {x.split(".")[0]: Chem.SDMolSupplier(os.path.join(self.args.inputs, x))
                              for x in self.input_files}

    @staticmethod
    def read_np_model(filename=None):
      """
      Reads and returns the scoring model, which has to be passed to the scoring functions.

      :param str filename: path to np scorer model
      :return:
      """
      if filename is None:
          filename = os.path.join(os.path.dirname(__file__), 'publicnp.model.gz')

      return pickle.load(gzip.open(filename))

    @staticmethod
    def read_fragment_scores(filename=None):
        """
        Reads and return the fragment scores which are used by the scorer

        :param filename:
        :return:
        """
        if filename is None:
            filename = os.path.join(os.path.dirname(__file__), 'fpscores.pkl.gz')
        _fscores = pickle.load(gzip.open(filename))

        return {j: i[0] for i in _fscores for j in i}

    def score_molecules(self):
        """
        Organises the iteration through compound sets

        :return:
        """
        dfs = {}

        for key, mols in self.molecule_sets.items():
            keys = []
            name = []
            smiles = []
            synthetic_accessibility = []
            natural_product_score = []
            print(key)
            for i, mol in enumerate(mols):
                if mol is None:
                    continue
                else:
                    keys.append(key)
                    name.append(i)
                    smiles.append(Chem.MolToSmiles(mol))
                    synthetic_accessibility.append(self.synthetic_accessibility.score(mol))
                    natural_product_score.append(self.natural_product_scorer.score(mol))

            df = pd.DataFrame({"Set": keys,
                               "Name": name,
                               "smiles": smiles,
                               "Synthetic Accessibility": synthetic_accessibility,
                               "Natural Product Likeness": natural_product_score})

            dfs[key] = df
        # df.to_csv(self.args.output_figure)
        return dfs

    def plot(self, dfs):
        sns.set_style("white")
        c = {"Abi": "Greens",
             "Drugs": "Blues",
             "NaturalProductsChEBI": "Reds"}

        for k, v in dfs.items():

            ax = sns.kdeplot(data=v['Synthetic Accessibility'], data2=v["Natural Product Likeness"], shade=False,
                             cmap=c[k])


            ax.set_xlim([0, 10])
            ax.set_ylim([-4, 4])

            plt.savefig(os.path.join(self.args.out_dir, k + ".png"))
            #plt.savefig(os.path.join(self.args.out_dir, "comparision.png"))
            plt.close()


if __name__ == '__main__':
    organiser = Organiser()
    dfs = organiser.score_molecules()
    organiser.plot(dfs)



