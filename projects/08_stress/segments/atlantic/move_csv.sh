#!/bin/bash

for file in $(ls *.txt)
do
  name=${file%"txt"}
  if [ ! -f ${name}"csv" ]; then
    echo ${name}
  fi
done 
