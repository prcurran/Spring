"""
Commandline program for Spring group to calculate the "Natural Product Likeness" scoring metric
"""

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import shutil


class Organiser(argparse.ArgumentParser):
    """
    class organising the natural product likeness calculation
    """

    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)

        self.add_argument(
            'inputs',
            help='Path to input directory. The directory should contain all molecule sets and no other data!'
        )

        self.add_argument(
            'output_figure',
            help='path to output figure'
        )

        self.add_argument(
            '-c', '--clean',
            default=True,
            help='If True the output .sdf files will be removed. (Recommended unless further analysis is required)'
        )

        self.args = self.parse_args()

        try:
            self.src = os.environ['NP_LIKENESS']

        except:
            raise EnvironmentError("Set np_likeness environment variable.  'export NP_LIKENESS=<path to NP-Likeness-2.1.jar>' ")

        self.data = {}
        self.outputs = {}
        self.supported_file_formats = [".sdf"]
        self.input_files = [x for x in os.listdir(self.args.inputs)
                            if os.path.splitext(x)[1] in self.supported_file_formats]

        self.output_directory = os.path.join(os.path.dirname(self.args.inputs), "output")

        if not os.path.exists(self.output_directory):
            os.mkdir(self.output_directory)

    def run_npl(self):
        """
        system call running the java file

        :return:
        """
        for i in self.input_files:
            input = os.path.join(self.args.inputs, i)
            output = os.path.join(self.output_directory, i)
            self.outputs.update({os.path.splitext(os.path.basename(i))[0]: output})

            cmd = "java -jar {} -in {} -out {} -intype sdf -outtype sdf".format(self.src, input, output)
            os.system(cmd)

    def get_scores(self):
        """
        extract scores from file

        :return:
        """
        for outname, out in self.outputs.items():

            f = open(out).readlines()
            scores = [float(f[i+1].strip('\n'))
                      for i, l in enumerate(f) if l == "> <NATURAL_PRODUCT_LIKENESS_SCORE>\n"
                      if len(f[i+1].strip('\n')) > 2]

            self.data.update({outname: pd.DataFrame({"{}".format(outname): scores})})

    def create_figure(self):
        """
        plot kde function for the score distributions

        :return:
        """
        for key, value in self.data.items():
            ax = sns.kdeplot(value[key])
        plt.savefig(self.args.output_figure)

    def cleanup(self):
        shutil.rmtree(self.output_directory)


def main():
    organiser = Organiser()
    organiser.run_npl()
    organiser.get_scores()
    organiser.create_figure()
    if organiser.args.clean is True:
        organiser.cleanup()


if __name__ == "__main__":
    main()