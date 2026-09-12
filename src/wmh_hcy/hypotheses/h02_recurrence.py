from ..analysis import run_survival


def run(data, cfg, folder):
    return run_survival(data, cfg, folder, "H2", risk=True)
