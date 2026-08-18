"""cdr_FS - feature selection by concentration/dose-response model fitting.

Selects features in high-content screening by fitting concentration/dose-response
models to earth mover's distance scores between control and treated cell populations.

The method was published in:

    Tome, M.; Jozef, B.; Mosimann, S. L.; Kosnik, M.; Schirmer, K.; Zupanic, A.
    A High-Content Imaging Pipeline to Investigate Subcytotoxic Effects in
    RTgill-W1 Cells. Environmental Science & Technology 2026, 60 (31), 21402-21416.
    https://doi.org/10.1021/acs.est.5c18316

Submodules import their scientific dependencies themselves; importing `cdr_fs` does
not pull in pandas, scipy or matplotlib, so the CLI stays usable for `--help` and for
configuration checks on a minimal install.
"""

__version__ = "0.1.0.dev0"

__all__ = ["__version__"]
