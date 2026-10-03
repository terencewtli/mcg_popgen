"""R1f (tables): genome-wide non-CpG P/D by substitution direction x element (gBGC control), 1 Mb block bootstrap.

Reads the per-chromosome R1f_gbgc_counts outputs. Compare S>W rows with the CpG-loss P/D of the same element in
R1e (CpG loss is S>W); W>S vs S>W P/D within an element measures gBGC there.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from R1e_pd_tables import boot_ci  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--counts-glob', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--n-boot', type=int, default=200)
    a = ap.parse_args()
    rng = np.random.default_rng(1)
    d = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(a.counts_glob))], ignore_index=True)
    d = d.groupby(['element', 'direction', 'block'])[['P', 'Ps', 'D']].sum().reset_index()
    rows = []
    for key, g in list(d.groupby(['element', 'direction'])) + [(('all', k), g) for k, g in d.groupby('direction')]:
        P, Ps, D = g.P.sum(), g.Ps.sum(), g.D.sum()
        lo, hi = boot_ci(g, a.n_boot, rng)
        rows.append(dict(element=key[0], direction=key[1], P=P, Ps=Ps, D=D, PD=P / D if D else np.nan,
                         PsD=Ps / D if D else np.nan, PD_lo=lo, PD_hi=hi))
    t = pd.DataFrame(rows)
    sw = t[t.direction == 'S>W'].set_index('element').PD
    ws = t[t.direction == 'W>S'].set_index('element').PD
    t['gbgc_SW_over_WS'] = t.element.map(sw / ws)
    bg = t[t.element == 'background'].set_index('direction').PD
    t['PD_rel_bg'] = t.PD / t.direction.map(bg)
    t.to_csv(a.out, sep='\t', index=False, float_format='%.4g')
    print(t.to_string(index=False, float_format=lambda v: f'{v:.3g}'), flush=True)


if __name__ == '__main__':
    main()
