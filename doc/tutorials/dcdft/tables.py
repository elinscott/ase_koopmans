from ase.collections import dcdft
import json
from pathlib import Path

data = json.loads(Path('fit.json').read_text())
for name in ['volume', 'B', 'Bp']:
    with open(name + '.csv', 'w') as f:
        print('# symbol, emt, exp, wien2k', file=f)
        for symbol, emt_dct in data.items():
            print(f'{symbol},', end=' ', file=f)
            dcdft_dct = dcdft.data[symbol]
            values = [emt_dct['emt_' + name],
                      dcdft_dct['exp_' + name],
                      dcdft_dct['wien2k_' + name]]
            print(', '.join(f'{value:.2f}'
                            for value in values), file=f)
