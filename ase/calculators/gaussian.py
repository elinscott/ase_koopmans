from shutil import which
from ase.io import read, write
from ase.calculators.calculator import FileIOCalculator, EnvironmentError


class Gaussian(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'dipole']
    command = 'GAUSSIAN < PREFIX.com > PREFIX.log'
    discard_results_on_any_change = True

    def __init__(self, *args, label='Gaussian', **kwargs):
        FileIOCalculator.__init__(self, *args, label='Gaussian', **kwargs)

    def calculate(self, *args, **kwargs):
        gaussians = ('g16', 'g09', 'g03')
        if 'GAUSSIAN' in self.command:
            for gau in gaussians:
                if which(gau):
                    self.command = self.command.replace('GAUSSIAN', gau)
                    break
            else:
                raise EnvironmentError('Missing Gaussian executable {}'
                                       .format(gaussians))

        FileIOCalculator.calculate(self, *args, **kwargs)

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write(self.label + '.com', atoms, properties=properties,
              format='gaussian-in', **self.parameters)

    def read_results(self):
        output = read(self.label + '.log', format='gaussian-out')
        self.calc = output.calc
        self.results = output.calc.results
