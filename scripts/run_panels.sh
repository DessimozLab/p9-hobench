#!/usr/bin/env bash
# Generic queue: MODE=full|subset PANELS="E P" scripts/run_panels.sh
# lane OF: orthofinder default -> of_flat ; lane FO: fastoma default (6 CPUs each lane).
# Progress: logs/wrapper_${MODE}_$(echo $PANELS|tr -d " ")_<lane>.log ; manifests in runs/<P>/full/<tool>_<variant>/.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
export PATH=$ROOT/conda_envs/ortho/bin:$PATH
export NXF_HOME=$ROOT/tools/nxf_home
PANELS=${PANELS:-"E P"}
MODE=${MODE:-subset}
mkdir -p logs

wait_panel() {  # tree + species.tsv + every proteome present, nothing modified for 3 min
  until python - "$1" "$MODE" <<'EOF2'
import sys, time
from pathlib import Path
from wrapper.common import load_panel
try:
    p = load_panel(sys.argv[1], sys.argv[2])
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
  if done_already "runs/$p/${MODE}/${tool}_${v}/manifest.json"; then echo "skip $tool $p $v (done)"; return 0; fi
  echo "[$(date -u +%FT%TZ)] START $tool $p $v"
  python -m wrapper.run --tool "$tool" --panel "$p" --mode ${MODE} --variant "$v" "$@"
  echo "[$(date -u +%FT%TZ)] END $tool $p $v rc=$?"
}

lane_of() {
  for p in $PANELS; do wait_panel "$p"; run orthofinder "$p" default --threads 6; done
  for p in $PANELS; do run of_flat "$p" default; done
}

lane_fo() {
  for p in $PANELS; do wait_panel "$p"; run fastoma "$p" default --threads 6 --mem-gb 30 --wait-hogmaps 180; done
}

lane_of > logs/wrapper_${MODE}_$(echo $PANELS|tr -d " ")_of.log 2>&1 &
lane_fo > logs/wrapper_${MODE}_$(echo $PANELS|tr -d " ")_fo.log 2>&1 &
wait
