from shutil import which
from ase.io import read, write
from ase.calculators.calculator import FileIOCalculator, EnvironmentError


class Gaussian(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'dipole']
    command = 'GAUSSIAN < PREFIX.com > PREFIX.log'

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
        write(self.label + '.com', atoms, properties, format='gaussian-in',
              **self.parameters)

    def read_output(self):
        output = read(self.label + '.log')
        self.calc = output.calc
        self.results = output.calc.results
