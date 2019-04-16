import os

from ccdc import io


base = "/home/pcurran/github_packages/np-likeness/data/input/3d"
fs = [os.path.join(base, f) for f in os.listdir(base) if f.endswith("sdf")]
mols = [io.MoleculeReader(fpath)[0] for fpath in fs]

with io.MoleculeWriter("abi_dos.sdf") as w:
    for m in mols:
        w.write(m)
