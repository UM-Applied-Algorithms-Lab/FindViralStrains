#!/bin/bash
Inputfile="$1"

# make a < .1 error file (1 path) (2nd column < .1)
awk '$2 < 0.1 {print > "1path.txt"} $2 >= 0.1 {print > "remaining.txt"}' $Inputfile
# make a < .1 error with remaining (2 path) and rest (3 path) (4th column < .1)
awk '$4 < 0.1 {print > "2path.txt"} $4 >= 0.1 {print > "3path.txt"}' remaining.txt

# sort each by their error
sort -k2,2 1path.txt -o 1path.txt
sort -k4,4 2path.txt -o 2path.txt
sort -k7,7 3path.txt -o 3path.txt

# put them together
cat 1path.txt 2path.txt 3path.txt > all.txt

# remove R2 lines
grep -v "_R2_" all.txt > allr1.txt
