
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter, _CSDDatabaseLocator
from ccdc.utilities import _private_importer
with _private_importer():
    import ChemicalAnalysisLib
    import ConformerGeneratorLib


def from_smiles(smiles, identifier=None, generate_initial_sites=True):
    """
    Create a :class:`ccdc.molecule.Molecule` from a SMILES string.
    *e.g.*::
         ethene = Molecule.from_smiles('C=C', 'Ethene')
    If ``identifier`` is not specified, the SMILES string will be used as the
    molecule identifier.
    :param smiles: str
    :param identifier: str
    :param generate_initial_sites: boolean - whether to include an initial guess at 3D coordinates
    :return: a :class:`ccdc.molecule.Molecule` instance with coordinates
    """
    try:
        if identifier is None:
            identifier = smiles

        if generate_initial_sites:
            parameter_files = _CSDDatabaseLocator.get_conformer_parameter_file_location()
            molmaker = ConformerGeneratorLib.MoleculeTo3D(parameter_files)
            mol = Molecule(identifier, molmaker.create_conformation(smiles))
        else:
            molmaker = ChemicalAnalysisLib.SMILESMoleculeMaker()
            mol = Molecule(identifier, _molecule=molmaker.siteless_atoms(smiles))
        return mol
    except RuntimeError:
        return None

def main():
    path = "/home/pcurran/github_packages/pyspring/np-likeness/scripts/NPs.smi"

    mols = [x.split("\t")[0] for x in open(path, "r").readlines()][1:]

    ccdc_mols = [from_smiles(x, identifier=str(i)) for i, x in enumerate(mols)]

    with MoleculeWriter("NP.sdf") as w:
        for m in [y for y in ccdc_mols if y is not None]:
            w.write(m)

    print len([y for y in ccdc_mols if y is not None])
    print len(mols)
if __name__ == "__main__":
    main()