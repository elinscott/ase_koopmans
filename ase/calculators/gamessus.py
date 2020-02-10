from ase.io import read, write
from ase.calculators.calculator import FileIOCalculator


class GAMESSUS(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'dipole']
    command = 'rungms PREFIX.inp > PREFIX.log'
    discard_results_on_any_change = True

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='gamessus', atoms=None, command=None, **kwargs):
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, command, **kwargs)

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write(self.label + '.inp', atoms, properties=properties,
              format='gamessus-in', **self.parameters)

    def read_results(self):
        output = read(self.label + '.log')
        self.calc = output.calc
        self.results = output.calc.results
