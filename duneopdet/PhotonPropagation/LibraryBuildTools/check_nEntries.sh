#!/bin/bash
where=$1
max=6000
if [ "$2" ]; then
  max=$2
fi
cr=0
for file in ${where}/*.root;
  do
    if [ $cr -lt $max ]; then
    toOpen=$file
    root -l -b -q "CountEntries.C+g(\"${toOpen}\")"
  fi
  cr=$((cr+1))
  done
