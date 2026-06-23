#!/usr/bin/env bash
# Filter pupmapper bedgraph files to only include positions with mappability > cutoff [default = 0.95]

BEDGRAPH_FILE=$1
CUTOFF=${2:-0.95}

awk -v cutoff="$CUTOFF" '$4 > cutoff' "$BEDGRAPH_FILE" > "${BEDGRAPH_FILE%.bedgraph}_filtered.${CUTOFF}.bedgraph"
