"""Python mirror of the E1 maps in scripts/slim/e1_passenger.slim (keep the two in lockstep; e1.py checks parity).

Sequence -> TF occupancy -> methylation -> CpG-transition gate. Fitness (selection models only) reads occupancy:
w = exp(-SEL * (1 - E)^2), E = mean occupancy of the two haplotypes.
"""
from dataclasses import dataclass, asdict

import numpy as np

MOTIF = np.array([3, 2, 0, 1, 3, 1, 0, 3])          # TGACTCAT, A=0 C=1 G=2 T=3; contains no CpG


@dataclass
class Params:
    N: int = 100
    L: int = 300
    ALPHA: float = 2e-5
    RCPG: float = 20.0
    MODE: int = 1            # 0 = context-only CpG hypermutation, 1 = methylation-gated
    SEL: float = 0.0         # 0 = neutral
    A0: float = 3.0
    A1: float = 0.0          # CpG-density feedback on methylation
    A2: float = 6.0          # TF protection
    W: int = 25
    LAMBDA: float = 25.0
    BETA: float = 2.0
    S0: float = 6.0
    BURN: int = 0
    MEAS: int = 0
    EVERY: int = 500
    NSAMP: int = 10

    @property
    def mstart(self):
        return self.L // 2 - 4

    @property
    def mcenter(self):
        return self.mstart + 3.5

    def slim_args(self):
        return asdict(self)


def occ_from_matches(s, P):
    return 1.0 / (1.0 + np.exp(-P.BETA * (s - P.S0)))


def occ(nuc, P):
    return occ_from_matches(np.sum(nuc[P.mstart:P.mstart + 8] == MOTIF), P)


def cpg_pos(nuc):
    return np.flatnonzero((nuc[:-1] == 1) & (nuc[1:] == 2))


def dist_to_motif(p, P):
    return np.maximum(0.0, np.abs(p - P.mcenter) - 4.0)


def meth_at(nuc, p, P, o=None):
    if len(p) == 0:
        return np.zeros(0)
    o = occ(nuc, P) if o is None else o
    rho = np.zeros(len(p))
    if P.A1 != 0:
        allp = cpg_pos(nuc)
        rho = (np.abs(allp[None, :] - p[:, None]) <= P.W).sum(1) - 1.0
    d = dist_to_motif(p, P)
    return 1.0 / (1.0 + np.exp(-(P.A0 - P.A1 * rho - P.A2 * o * np.exp(-d / P.LAMBDA))))


def gate(m, P):
    """Acceptance probability of a proposed CpG transition: rate ALPHA * (1 + (RCPG - 1) m) out of RCPG * ALPHA."""
    return (1.0 + (P.RCPG - 1.0) * m) / P.RCPG


def sel_for_two_ns(two_ns_target, N, P):
    """SEL such that a heterozygote carrying one motif mismatch (8/8 + 7/8) has 2N*s = two_ns_target vs 8/8 + 8/8."""
    e0 = occ_from_matches(8, P)
    e1 = 0.5 * (occ_from_matches(8, P) + occ_from_matches(7, P))
    dq = (1 - e1) ** 2 - (1 - e0) ** 2
    s = two_ns_target / (2 * N)
    return -np.log(1 - s) / dq


def selection_table(P):
    """2N*s (relative to the 8/8 homozygote) for the genotypes one or two mismatches away."""
    o = {k: occ_from_matches(k, P) for k in range(9)}
    w = lambda a, b: np.exp(-P.SEL * (1 - 0.5 * (o[a] + o[b])) ** 2)
    rows = {'8/7': w(8, 7), '7/7': w(7, 7), '8/6': w(8, 6), '8/5': w(8, 5), '8/0': w(8, 0), '0/0': w(0, 0)}
    return {k: 2 * P.N * (1 - v / w(8, 8)) for k, v in rows.items()}


# ---------------------------------------------------------------------------------------------------------------
# E2 readouts (scripts/slim/e2_stab.slim): GMODE 1 = P, g = O; GMODE 2 = C, g = exp(-R / R0), R = sum of m over the
# locus. Fitness w = exp(-KAPPA * (E - ESTAR)^2), E = mean readout of the two haplotypes.

def meth_load(nuc, P):
    return float(np.sum(meth_at(nuc, cpg_pos(nuc), P)))


def readout(nuc, P, gmode, R0=3.0):
    if gmode == 1:
        return float(occ(nuc, P))
    return float(np.exp(-meth_load(nuc, P) / R0))


# ---------------------------------------------------------------------------------------------------------------
# E3 (scripts/slim/e3_combined.slim): motif choice and combined readouts.
MOTIFS = {0: MOTIF, 1: np.array([3, 2, 0, 1, 2, 3, 1, 0])}     # 0 TGACTCAT; 1 TGACGTCA (CRE, CpG at motif pos 3-4)


def occ_m(nuc, P, motif):
    """CRE: the motif CpG (MSTART+3, +4) is required for binding; O = 0 without it (mirrors e3_combined.slim)."""
    if motif is MOTIFS[1] and not (nuc[P.mstart + 3] == 1 and nuc[P.mstart + 4] == 2):
        return 0.0
    return float(occ_from_matches(np.sum(nuc[P.mstart:P.mstart + 8] == motif), P))


def meth_at_m(nuc, p, P, motif):
    if len(p) == 0:
        return np.zeros(0)
    o = occ_m(nuc, P, motif)
    rho = np.zeros(len(p))
    if P.A1 != 0:
        allp = cpg_pos(nuc)
        rho = (np.abs(allp[None, :] - p[:, None]) <= P.W).sum(1) - 1.0
    d = dist_to_motif(p, P)
    return 1.0 / (1.0 + np.exp(-(P.A0 - P.A1 * rho - P.A2 * o * np.exp(-d / P.LAMBDA))))


def readout_e3(nuc, P, gmode, motif_id, R0=3.0):
    """GMODE 1: O; 2: exp(-R/R0); 3: O exp(-R/R0); 4: O (1 - m_c), m_c = methylation of the CpG at MSTART+3."""
    motif = MOTIFS[motif_id]
    o = occ_m(nuc, P, motif)
    if gmode == 1:
        return o
    if gmode == 4:
        c = P.mstart + 3
        mc = float(meth_at_m(nuc, np.array([c]), P, motif)[0]) if (nuc[c] == 1 and nuc[c + 1] == 2) else 0.0
        return o * (1.0 - mc)
    R = float(np.sum(meth_at_m(nuc, cpg_pos(nuc), P, motif)))
    return o * np.exp(-R / R0) if gmode == 3 else float(np.exp(-R / R0))
