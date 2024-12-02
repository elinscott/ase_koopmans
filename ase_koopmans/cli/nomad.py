import os
import os.path as op
import subprocess


class CLICommand:
    """Upload files to NOMAD.

    Upload all data within specified folders to the Nomad repository
    using authentication token given by the --token option or,
    if no token is given, the token stored in ~/.ase_koopmans/nomad-token.

    To get an authentication token, you create a Nomad repository account
    and use the 'Uploads' button on that page while logged in:

      https://repository.nomad-coe.eu/
    """

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('folders', nargs='*', metavar='folder')
        parser.add_argument('-t', '--token',
                            help='Use given authentication token and save '
                            'it to ~/.ase_koopmans/nomad-token unless '
                            '--no-save-token')
        parser.add_argument('-n', '--no-save-token', action='store_true',
                            help='do not save the token if given')
        parser.add_argument('-0', '--dry-run', action='store_true',
                            help='print command that would upload files '
                            'without uploading anything')

    @staticmethod
    def run(args):
        dotase_koopmans = op.expanduser('~/.ase_koopmans')
        tokenfile = op.join(dotase_koopmans, 'nomad-token')

        if args.token:
            token = args.token
            if not args.no_save_token:
                if not op.isdir(dotase_koopmans):
                    os.mkdir(dotase_koopmans)
                with open(tokenfile, 'w') as fd:
                    print(token, file=fd)
                os.chmod(tokenfile, 0o600)
                print('Wrote token to', tokenfile)
        else:
            try:
                with open(tokenfile) as fd:
                    token = fd.readline().strip()
            except OSError as err:  # py2/3 discrepancy
                from ase_koopmans.cli.main import CLIError
                msg = ('Could not find authentication token in {}.  '
                       'Use the --token option to specify a token.  '
                       'Original error: {}'
                       .format(tokenfile, err))
                raise CLIError(msg)

        cmd = ('tar cf - {} | '
               'curl -XPUT -# -HX-Token:{} '
               '-N -F file=@- http://nomad-repository.eu:8000 | '
               'xargs echo').format(' '.join(args.folders), token)

        if not args.folders:
            print('No folders specified -- another job well done!')
        elif args.dry_run:
            print(cmd)
        else:
            print('Uploading {} folder{} ...'
                  .format(len(args.folders),
                          's' if len(args.folders) != 1 else ''))
            subprocess.check_call(cmd, shell=True)
