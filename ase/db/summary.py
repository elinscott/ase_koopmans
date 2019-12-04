
from ase.geometry import cell_to_cellpar

if 0:
    def write(self):
        print(self.formula + ':')
        for headline, columns in self.layout:
            blocks = columns[0]
            if len(columns) == 2:
                blocks += columns[1]
            print((' ' + headline + ' ').center(78, '='))
            for block in blocks:
                if block['type'] == 'table':
                    rows = block['rows']
                    if not rows:
                        print()
                        continue
                    rows = [block['header']] + rows
                    widths = [max(len(row[n]) for row in rows)
                              for n in range(len(rows[0]))]
                    for row in rows:
                        print('|'.join('{:{}}'.format(word, width)
                                       for word, width in zip(row, widths)))
                    print()
                elif block['type'] == 'figure':
                    print(block['filename'])
                    print()
                elif block['type'] == 'cell':
                    print('Unit cell in Ang:')
                    print('axis|periodic|          x|          y|          z')
                    c = 1
                    fmt = '   {0}|     {1}|{2[0]:>11}|{2[1]:>11}|{2[2]:>11}'
                    for p, axis in zip(self.row.pbc, self.cell):
                        print(fmt.format(c, [' no', 'yes'][p], axis))
                        c += 1
                    print('Lengths: {:>10}{:>10}{:>10}'
                          .format(*self.lengths))
                    print('Angles:  {:>10}{:>10}{:>10}\n'
                          .format(*self.angles))

        if self.stress:
            print('Stress tensor (xx, yy, zz, zy, zx, yx) in eV/Ang^3:')
            print('   ', self.stress, '\n')

        if self.dipole:
            print('Dipole moment in e*Ang: ({})\n'.format(self.dipole))

        if self.constraints:
            print('Constraints:', self.constraints, '\n')

        if self.data:
            print('Data:', self.data, '\n')
