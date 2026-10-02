#!/bin/bash
# Sync curated docs/code/small results from the working directory into the git mirror (github/mcg_popgen/), then
# optionally commit. Content flows working dir -> mirror only; never hand-edit mirror copies of synced files.
# Exception: md/JOURNAL.md, md/RESULTS.md, md/PROGRESS.md live ONLY in the mirror (edit them there).
# realdata/ is owned by the cluster session and is never touched here (no rsync targets it).
#
# Usage: bash scripts/sync_to_github.sh ["commit message"]
# With no message, stages and shows the status but does not commit. Never pushes.
#
# Never mirrored: raw simulation output (sim/, *.trees, *.vcf, *.npz) — regenerate from scripts + seeds.
# Small tables go through csv/ or tsv/; figures through pdf/ or png/ (any single file > 20 MB is skipped).

set -euo pipefail

PROJDIR=/Users/terenceli/claude/project_ideas/mcg_popgen
GH=$PROJDIR/github/mcg_popgen
MSG=${1:-}

[ -d "$GH/.git" ] || { echo "no git repo at $GH" >&2; exit 1; }
# the cluster session commits realdata/ to the same repo (md/ENDPOINT.md); pick its commits up first
git -C "$GH" pull --ff-only --quiet || { echo "git pull --ff-only failed; resolve in $GH first" >&2; exit 1; }

[ -f "$PROJDIR/README.md" ] && cp "$PROJDIR/README.md" "$GH/README.md"
mkdir -p "$GH/md" "$GH/scripts"

# notes (md/): everything except the mirror-only files (excluded files are never deleted by --delete)
rsync -a --delete --exclude='JOURNAL.md' --exclude='RESULTS.md' --exclude='PROGRESS.md' \
    --include='*/' --include='*.md' --exclude='*' "$PROJDIR/md/" "$GH/md/"

# code: SLiM/Eidos, Python, shell, R, yaml, md — every scripts/ subdirectory
rsync -a --delete --include='*/' --include='*.slim' --include='*.eidos' --include='*.py' --include='*.sh' \
    --include='*.R' --include='*.yaml' --include='*.md' \
    --exclude='__pycache__/' --exclude='*' --prune-empty-dirs "$PROJDIR/scripts/" "$GH/scripts/"
if [ -d "$PROJDIR/tests" ]; then
    rsync -a --delete --include='*/' --include='*.py' --include='*.slim' --exclude='__pycache__/' --exclude='*' \
        --prune-empty-dirs "$PROJDIR/tests/" "$GH/tests/"
fi

# LaTeX sources (compiled PDFs go through pdf/)
if [ -d "$PROJDIR/tex" ]; then
    rsync -a --delete --include='*/' --include='*.tex' --include='*.bib' --exclude='*' --prune-empty-dirs \
        "$PROJDIR/tex/" "$GH/tex/"
fi

# small tables (summary statistics, benchmarks)
for d in csv tsv; do
    [ -d "$PROJDIR/$d" ] || continue
    rsync -a --delete --include='*/' --include='*.tsv' --include='*.csv' --exclude='*' --max-size=20m \
        --prune-empty-dirs "$PROJDIR/$d/" "$GH/$d/"
done

# figures
for d in pdf png; do
    [ -d "$PROJDIR/$d" ] || continue
    rsync -a --delete --include='*/' --include="*.$d" --exclude='*' --max-size=20m --prune-empty-dirs \
        "$PROJDIR/$d/" "$GH/$d/"
done

# notebooks (checkpoints excluded); outputs kept — sims are small, so notebooks should be too
if [ -d "$PROJDIR/ipynb" ]; then
    rsync -a --delete --include='*/' --include='*.ipynb' --exclude='.ipynb_checkpoints/' --exclude='*' \
        --max-size=20m --prune-empty-dirs "$PROJDIR/ipynb/" "$GH/notebooks/"
fi

cd "$GH"
git add -A
git status --short
if [ -n "$MSG" ]; then
    git commit -m "$MSG"
else
    echo "(no commit message given — staged only)"
fi
