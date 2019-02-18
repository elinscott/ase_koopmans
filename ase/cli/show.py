from ase.io.jsonio import read_json


class CLICommand:
    """Describe and/or visualize ASE data.

    Examples:

        $ # Show band path stored in bandpath.json
        $ ase show bandpath.json
        $ # Show band structure stored in bs.json
        $ ase show bs.json
      """

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('file', nargs='+')

    @staticmethod
    def run(args, parser):
        for fname in args.file:
            obj = read_json(fname)
            objtype = obj.ase_objtype
            print(objtype)
            print(obj)

            # plot() should be uniform among plottable objects
            obj.plot()

        import matplotlib.pyplot as plt
        plt.show()
