"""Quantum ESPRESSO Calculator

export ASE_ESPRESSO_KCP_COMMAND="/path/to/kcp.x -in PREFIX.cpi > PREFIX.cpo"

N.B. the extensions must be .cpi and .cpo

Run kcp.x jobs.
"""


from ._espresso import EspressoParent


class Espresso_kcp(EspressoParent):
    ext_in = '.cpi'
    ext_out = '.cpo'
    implemented_properties = ['energy']

    # Default command does not use parallelism and assumes kcp.x is on the user's path
    command = 'kcp.x -in PREFIX.cpi > PREFIX.cpo 2>&1'

    def __init__(self, *args, **kwargs):
        kwargs['label'] = 'kc_ham'
        super().__init__(*args, **kwargs)
