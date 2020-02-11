import warnings
import os
from subprocess import Popen

from ase.io import read, write
from ase.calculators.calculator import FileIOCalculator


class GAMESSUS(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'dipole']
    command = 'rungms PREFIX.inp > PREFIX.log'
    discard_results_on_any_change = True

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='gamessus', atoms=None, command=None, userscr=None,
                 **kwargs):
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, command, **kwargs)
        self.userscr = userscr

    def calculate(self, *args, **kwargs):
        if self.userscr is None:
            if 'rungms' in self.command:
                self.userscr = self.get_userscr()

        if self.userscr is None:
            warnings.warn("Could not determine USERSCR! "
                          "GAMESS may refuse to run more than once for "
                          "this job. Please pass userscr to the GAMESSUS "
                          "Calculator if you run into problems!")
        else:
            self.clean_userscr()

        FileIOCalculator.calculate(self, *args, **kwargs)

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write(self.label + '.inp', atoms, properties=properties,
              format='gamessus-in', **self.parameters)

    def read_results(self):
        output = read(self.label + '.log')
        self.calc = output.calc
        self.results = output.calc.results

    def get_userscr(self):
        """Runs rungms without an input file to determine USERSCR"""
        # I'm probably overthinking this, but what if prefix_test.inp exists?
        # Keep appending _test until we get a name that doesn't exist.
        prefix = self.prefix + '_test'
        while os.path.exists(prefix + '.inp'):
            prefix += '_test'

        command = self.command.replace('PREFIX', prefix)
        proc = Popen(command, shell=True, cwd=self.directory)
        err = proc.wait()

        with open(prefix + '.log') as f:
            for line in f:
                if line.startswith('GAMESS supplementary output files'):
                    # in case USERSCR has spaces in it...
                    return ' '.join(line.split(' ')[8:]).strip()
        return None

    def clean_userscr(self):
        """Backs up any conflicting files in USERSCR"""
        assert self.userscr is not None
        for fname in os.listdir(self.userscr):
            tokens = fname.split('.')
            if tokens[0] == self.prefix and tokens[-1] != 'bak':
                fold = os.path.join(self.userscr, fname)
                os.rename(fold, fold + '.bak')
