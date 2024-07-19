#! /bin/bash

for data in ind/ anti/ corr/
do
    ./main PDBB $data 10000 6 0.5 6 >> res.txt
    for theta in 10 20 40 80 160 320 640 1280
    do
        ./main PRE $data 10000 6 0.5 6 $theta >> res.txt
    done
done