"""Single-lineage oracle for the E0 neutral nucleotide model.

Under neutrality, a haplosome sampled at generation t is the ancestral sequence evolved for t generations along one
lineage under the mutation process alone (drift reshuffles which lineage is sampled, not what a lineage experiences).
So E[CpG count at t] in SLiM must equal E[CpG count] of a continuous-time Markov chain on sequences, started from the
same ancestral sequence, with SLiM's per-generation rates. Discrete vs continuous time differs at O(rate^2).

The rate model mirrors scripts/slim/e0_neutral.slim exactly, including SLiM's convention that bases off either
chromosome end count as "A" for context. CpG transitions (C->T at a CpG C, G->A at a CpG G) have rate
gate * rcpg * alpha; every other change has rate alpha.
"""
import numpy as np

NUC = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def rate_matrix(alpha, rcpg, gate=1.0):
    """64x4 array: rate of the central base of trinucleotide (a, b, c), indexed a*16 + b*4 + c, mutating to each base."""
    mm = np.full((64, 4), alpha)
    for i in range(64):
        mm[i, (i // 4) % 4] = 0.0
    for a in range(4):
        mm[a * 16 + 1 * 4 + 2, 3] = gate * rcpg * alpha   # aCG: C -> T
        mm[1 * 16 + 2 * 4 + a, 0] = gate * rcpg * alpha   # CGa: G -> A
    return mm


def cpg_count(seq):
    return int(np.sum((seq[:-1] == 1) & (seq[1:] == 2)))


def simulate_lineage(anc, mm, times, rng):
    """Gillespie along one lineage; returns the CpG count at each time in `times` (sorted, in generations)."""
    s = anc.copy()
    L = len(s)
    pad = np.concatenate([[0], s, [0]])                     # off-end bases are "A"

    def site_rates(idx):
        tri = pad[idx] * 16 + pad[idx + 1] * 4 + pad[idx + 2]
        return mm[tri]                                      # (len(idx), 4)

    R = site_rates(np.arange(L))                           # per-site, per-target rates
    tot = R.sum(1)
    out = np.empty(len(times), dtype=int)
    t, k = 0.0, 0
    while k < len(times):
        T = tot.sum()
        t += rng.exponential(1.0 / T)
        while k < len(times) and t > times[k]:
            out[k] = cpg_count(s)
            k += 1
        if k == len(times):
            break
        i = np.searchsorted(np.cumsum(tot), rng.random() * T)
        i = min(i, L - 1)
        new = rng.choice(4, p=R[i] / tot[i])
        s[i] = new
        pad[i + 1] = new
        nb = np.arange(max(0, i - 1), min(L, i + 2))
        R[nb] = site_rates(nb)
        tot[nb] = R[nb].sum(1)
    return out


def expected_cpg_trajectory(anc_str, alpha, rcpg, gate, times, n_lineages, seed):
    anc = np.array([NUC[c] for c in anc_str.strip()])
    mm = rate_matrix(alpha, rcpg, gate)
    rng = np.random.default_rng(seed)
    traj = np.array([simulate_lineage(anc, mm, np.asarray(times, float), rng) for _ in range(n_lineages)])
    return traj.mean(0), traj.std(0, ddof=1) / np.sqrt(n_lineages)


# ---------------------------------------------------------------------------------------------------------------
# E1 oracle: one lineage under the full sequence -> occupancy -> methylation -> gated-mutation map (model.py).
# Exact for neutral models (SEL = 0): each haplotype's mutation process depends only on its own sequence.
# With freeze_motif=True the motif bases never mutate, approximating strong selection on occupancy (P models);
# SLiM - oracle then measures the effect of motif polymorphism / mismatch load.

def _e1_rates(s, P, mm, freeze_motif):
    import model
    pad = np.concatenate([[0], s, [0]])
    tri = pad[:-2] * 16 + pad[1:-1] * 4 + pad[2:]
    R = mm[tri].copy()
    if P.MODE == 1:
        p = model.cpg_pos(s)
        if len(p):
            g = model.gate(model.meth_at(s, p, P), P)
            R[p, 3] *= g          # C -> T at the CpG C
            R[p + 1, 0] *= g      # G -> A at the CpG G
    if freeze_motif:
        R[P.mstart:P.mstart + 8] = 0.0
    return R


def e1_lineage(anc, P, sample_times, rng, freeze_motif=False):
    """Returns per-position sums over the sample times: cpg, C, G, msum, plus summed occupancy and #samples."""
    import model
    mm = rate_matrix(P.ALPHA, P.RCPG)
    s = anc.copy()
    L = len(s)
    acc = dict(cpg=np.zeros(L), c=np.zeros(L), g=np.zeros(L), msum=np.zeros(L), osum=0.0, n=0)
    R = _e1_rates(s, P, mm, freeze_motif)
    tot = R.sum(1)
    t, k = 0.0, 0
    while k < len(sample_times):
        T = tot.sum()
        t += rng.exponential(1.0 / T)
        while k < len(sample_times) and t > sample_times[k]:
            p = model.cpg_pos(s)
            acc['cpg'][p] += 1
            acc['c'][s == 1] += 1
            acc['g'][s == 2] += 1
            acc['msum'][p] += model.meth_at(s, p, P)
            acc['osum'] += model.occ(s, P)
            acc['n'] += 1
            k += 1
        if k == len(sample_times):
            break
        i = min(np.searchsorted(np.cumsum(tot), rng.random() * T), L - 1)
        s[i] = rng.choice(4, p=R[i] / tot[i])
        R = _e1_rates(s, P, mm, freeze_motif)    # full recompute: occupancy and CpG density are non-local
        tot = R.sum(1)
    return acc
