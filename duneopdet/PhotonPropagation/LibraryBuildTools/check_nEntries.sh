#!/bin/bash
where=$1
cr=0
for file in ${where}/*.root;
  do
    if [ $cr -lt 3 ]; then
    toOpen=$file
    root -b -q "CountEntries.C(\"${toOpen}\")"
  fi
  cr=$((cr+1))
  done
