#!/bin/bash

where=$1
BatchID=$2
echo "Assembling library for files from $where into $2"

list=temp${BatchID}_FileList.txt

#Produce the list of files to combine
ls -1 $where/root/*.root > $list

root -b -q "AssembleSingleFile.C+g(\"$list\", \"\",\"lib_$BatchID.root\")"

cp lib_$BatchID.root $where

#clean up temporaryfiles
rm $list

