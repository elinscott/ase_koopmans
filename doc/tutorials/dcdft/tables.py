import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from ase.collections import dcdft
from ase.eos import birchmurnaghan
from ase.units import kJ
from ase.utils.dcdft import delta

# Read EMT data:
data = json.loads(Path('fit.json').read_text())
# Insert values from experiment and WIEN2k:
for symbol in data:
    dcdft_dct = dcdft.data[symbol]
    dcdft_dct['exp_B'] *= 1e-24 * kJ
    dcdft_dct['wien2k_B'] *= 1e-24 * kJ
    data[symbol].update(dcdft_dct)

for name in ['volume', 'B', 'Bp']:
    with open(name + '.csv', 'w') as f:
        print('# symbol, emt, exp, wien2k', file=f)
        for symbol, dct in data.items():
            values = [dct[code + '_' + name]
                      for code in ['emt', 'exp', 'wien2k']]
            if name == 'B':
                values = [val * 1e24 / kJ for val in values]
            print(f'{symbol},',
                  ', '.join(f'{value:.2f}' for value in values),
                  file=f)

with open('delta.csv', 'w') as f:
    print('# symbol, emt-exp, emt-wien2k, exp-wien2k', file=f)
    for symbol, dct in data.items():
        # Get v0, B, Bp:
        emt, exp, wien2k = [(dct[code + '_volume'],
                             dct[code + '_B'],
                             dct[code + '_Bp'])
                            for code in ['emt', 'exp', 'wien2k']]
        print(f'{symbol},',
              '{:.1f}, {:.1f}, {:.1f}'.format(delta(*emt, *exp) * 1000,
                                              delta(*emt, *wien2k) * 1000,
                                              delta(*exp, *wien2k) * 1000),
              file=f)

        if symbol == 'Pt':
            for (v0, B, Bp), code in [(emt, 'EMT'),
                                      (exp, 'experiment'),
                                      (wien2k, 'WIEN2k')]:
                v = np.linspace(0.94 * v0, 1.06 * v0)
                plt.plot(v, birchmurnaghan(v, 0.0, B, Bp, v0), label=code)
            plt.legend()
            plt.xlabel('volume [Ang^3]')
            plt.ylabel('energy [eV/atom]')
            plt.savefig('Pt.png')
