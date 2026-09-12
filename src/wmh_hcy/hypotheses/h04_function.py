from ..analysis import run_functional


def run(data, cfg, folder, t1=False):
    return run_functional(data, cfg, folder, "functional_t1" if t1 else "functional")
