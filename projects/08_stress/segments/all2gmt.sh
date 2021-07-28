#!/bin/bash

for file in $(ls *.kml)
do
  gmt kml2gmt $file -V > ${file%.kml}".txt"
done
