#!/bin/bash

# Gokberk Alagoz
# 10.05.2021

# Sort bed files based on chr position and startPos
echo "Enter the path to .beds: "
read folder

for f in $folder/*.bed; do
	sort -V -k1,1 -k2,2 $f > ${f%.bed}.sorted.bed
done