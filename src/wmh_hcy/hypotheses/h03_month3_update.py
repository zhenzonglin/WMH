from ..analysis import run_survival
from ..design import Design


def run(data, cfg, folder):
    return run_survival(data, cfg, folder, "H3", spec=Design(month3=True), kind="month3", risk=True)
