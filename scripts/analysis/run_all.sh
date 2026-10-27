#!/bin/bash
# Regenerate results/analysis/ from results/<P>_<mode>/ and runs/ (read-only). Run from the repo root.
set -euo pipefail
cd "$(dirname "$0")"
for pm in Y:full K:full E:subset P:subset V:subset; do
  P=${pm%:*}; M=${pm#*:}
  [ -f ../../results/${P}_${M}/per_family.tsv ] && python extract_nodes.py --panel $P --mode $M
done
python main_table.py; python depth_profile.py; python phenomena.py; python difficulty.py; python error_anatomy.py; python robustness.py
python ../confirmatory.py --spec Y:full K:full E:subset P:subset V:subset --out ../../results/analysis/confirmatory_current.tsv
