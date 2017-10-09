#!/bin/bash
where=$1
max=$2
cr=0
for file in ${where}/*.root;
  do
    if [ $cr -lt $max ]; then
    toOpen=$file
    root -l -b -q "CountEntries.C+g(\"${toOpen}\")"
  fi
  cr=$((cr+1))
  done
