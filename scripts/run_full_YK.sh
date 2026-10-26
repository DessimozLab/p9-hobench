#!/usr/bin/env bash
# Full-mode tool runs for panels Y and K (DESIGN §7 E1, SYNTHESIS E4/E7).
# Two lanes in parallel (12 cores, ~50 GB RAM):
#   lane OF : orthofinder default Y, K (6 threads) -> of_flat Y, K -> orthofinder dendroblast Y, K
#   lane FO : fastoma Y, K (6 CPUs, 30 GB cap; waits up to 3 h for precomputed hogmaps)
# Progress: logs/wrapper_YK_<lane>.log ; manifests in runs/<P>/full/<tool>_<variant>/.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
export PATH=$ROOT/conda_envs/ortho/bin:$PATH
export NXF_HOME=$ROOT/tools/nxf_home
PANELS=${PANELS:-"Y K"}
mkdir -p logs

wait_panel() {  # tree + species.tsv + every proteome present, nothing modified for 3 min
  until python - "$1" <<'EOF2'
import sys, time
from pathlib import Path
from wrapper.common import load_panel
try:
    p = load_panel(sys.argv[1], "full")
except Exception:
    sys.exit(1)
files = [p.tree_path, p.dir / "species.tsv"] + list(p.fastas.values())
sys.exit(0 if time.time() - max(f.stat().st_mtime for f in files) > 180 else 1)
EOF2
  do sleep 60; done
  echo "[$(date -u +%FT%TZ)] panel $1 ready"
}

done_already() { python - "$1" <<'EOF'
import json,sys
try: sys.exit(0 if json.load(open(sys.argv[1]))["status"]=="done" else 1)
except Exception: sys.exit(1)
EOF
}

run() {  # tool panel variant [extra args]
  local tool=$1 p=$2 v=$3; shift 3
  if done_already "runs/$p/full/${tool}_${v}/manifest.json"; then echo "skip $tool $p $v (done)"; return 0; fi
  echo "[$(date -u +%FT%TZ)] START $tool $p $v"
  python -m wrapper.run --tool "$tool" --panel "$p" --mode full --variant "$v" "$@"
  echo "[$(date -u +%FT%TZ)] END $tool $p $v rc=$?"
}

lane_of() {
  for p in $PANELS; do wait_panel "$p"; run orthofinder "$p" default --threads 6; done
  for p in $PANELS; do run of_flat "$p" default; done
  for p in $PANELS; do run orthofinder "$p" dendroblast --threads 6; run of_flat "$p" dendroblast; done
}

lane_fo() {
  for p in $PANELS; do wait_panel "$p"; run fastoma "$p" default --threads 6 --mem-gb 30 --wait-hogmaps 180; done
}

lane_of > logs/wrapper_YK_of.log 2>&1 &
lane_fo > logs/wrapper_YK_fo.log 2>&1 &
wait
