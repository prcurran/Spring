import os

mols = "/home/pcurran/github_packages/Spring/data/input"
out = "/home/pcurran/github_packages/Spring/data/output"
cmd = "python /home/pcurran/github_packages/Spring/src/sa_np_scorer.py {} {}".format(mols, out)

os.system(cmd)