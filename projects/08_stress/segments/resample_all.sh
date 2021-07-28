#!/bin/bash

for file in $(ls *.txt)
do
  gmt sample1d $file -T1m -AR > temp
  mv temp $file
done
