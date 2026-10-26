#!/usr/bin/env bash
# List files that appeared in requested/ besides REQUESTED.md, with size and sha256 (run at each milestone).
cd "$(dirname "$0")/../requested" || exit 1
find . -maxdepth 2 -type f ! -name REQUESTED.md ! -name RECEIVED.tsv -printf '%TY-%Tm-%Td %TH:%TM\t%s\t%p\n' | sort
