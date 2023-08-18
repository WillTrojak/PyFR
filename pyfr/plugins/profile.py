from pyfr.plugins.base import BaseSolnPlugin


class ProfilePlugin(BaseSolnPlugin):
    name = 'profile'
    systems = ['*']
    formulations = ['dual', 'std']

    def __init__(self, intg, cfgsect, prefix):
        super().__init__(intg, cfgsect, prefix)

        self.n_start = self.cfg.getint(cfgsect, 'n-start', 1)
        self.n_stop = self.cfg.getint(cfgsect, 'n-stop', -1)

    def __call__(self, intg):
        if intg.nacptsteps == self.n_start:
            intg.system.backend.start_profile_trace()

        if intg.nacptsteps == self.n_stop:
            intg.system.backend.stop_profile_trace()
