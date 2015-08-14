#!/bin/bash

where=$1
BatchID=$2
echo "Assembling library for files from $where into $2"

list=temp${BatchID}_FileList.txt

#Produce the list of files to combine
ls $where/root/*.root > $list

#Make a string to feed to roor
#echo ".L AssembleSingleFile.C" >> rootstring.txt
#echo "AssembleSingleFile(\"temp${BatchID}_FileList.txt\", \"root/Files$BatchID/\",\"photlibrary/lib$BatchID.root\")" >> rootstring.txt

#Make root execute this command
#cat rootstring.txt | root -l

root -b -q "AssembleSingleFile.C(\"$list\", \"\",\"$where/lib_$BatchID.root\")"

#clean up temporaryfiles
rm $list
#rm rootstring.txt

