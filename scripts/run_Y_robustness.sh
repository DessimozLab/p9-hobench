#!/usr/bin/env bash
# E2 robustness on Y full (seed 20261026): waits for the E/P subset queue, then OF+FastOMA on
# Y_alttree, Y_frag, Y_contam, Y_taxrm (priority d,a,c,e), then OF clean replicate (rep2) on Y.
cd "$(dirname "$0")/.."
[ -n "${1:-}" ] && while kill -0 "$1" 2>/dev/null; do sleep 60; done  # optional: PID of a queue to wait for
MODE=full PANELS="Y_alttree Y_frag Y_contam Y_taxrm" WAIT_HOGMAPS=0 scripts/run_panels.sh
export PATH=$PWD/conda_envs/ortho/bin:$PATH
python -m wrapper.run --tool orthofinder --panel Y --mode full --variant rep2 --threads 6 > logs/wrapper_Y_rep2.log 2>&1
