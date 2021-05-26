#!/usr/bin/env bash

ssp=$1
grid=$2

gmt grdtrack ${ssp} -G${grid} -C20k/0.5k/1k+v -Sm+sSTACKS/${ssp%.xy}_stack.txt > STACKS/table.txt


## OUTPUT
## The leading column holds cross distance, while the first four columns in a
## group hold stacked value, deviation, min value, and max value, respectively.
## If method is one of a|m|p then we also write the lower and upper confidence 
## bounds (see +c).
