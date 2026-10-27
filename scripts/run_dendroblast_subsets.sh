#!/usr/bin/env bash
# OrthoFinder -M dendroblast + of_flat on E, P, V subsets (third method). Optional arg: PID to wait for.
cd "$(dirname "$0")/.."
[ -n "${1:-}" ] && while kill -0 "$1" 2>/dev/null; do sleep 60; done
export PATH=$PWD/conda_envs/ortho/bin:$PATH
for p in E P V; do
  for t in orthofinder of_flat; do
    echo "[$(date -u +%FT%TZ)] START $t $p dendroblast"
    python -m wrapper.run --tool $t --panel $p --mode subset --variant dendroblast --threads 12
    echo "[$(date -u +%FT%TZ)] END $t $p dendroblast rc=$?"
  done
done
