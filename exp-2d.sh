#! /bin/bash

for data in ind/ anti/ corr/
do
    for n in 1000 10000 50000 100000 500000 1000000
    do
        timeout 1h ./main BSL $data $n 2 0.5 6 >> res.txt
        timeout 1h ./main DP $data $n 2 0.5 6 >> res.txt
        timeout 1h ./main DP+ $data $n 2 0.5 6 >> res.txt
        timeout 1h ./main GDY $data $n 2 0.5 6 0.01 >> res.txt
        timeout 1h ./main PDBB $data $n 2 0.5 6 >> res.txt
        timeout 1h ./main PRE $data $n 2 0.5 6 160 >> res.txt
    done

    for p in 0.2 0.35 0.5 0.65 0.8
    do
        timeout 1h ./main BSL $data 10000 2 $p 6 >> res.txt
        timeout 1h ./main DP $data 10000 2 $p 6 >> res.txt
        timeout 1h ./main DP+ $data 10000 2 $p 6 >> res.txt
        timeout 1h ./main GDY $data 10000 2 $p 6 0.01 >> res.txt
        timeout 1h ./main PDBB $data 10000 2 $p 6 >> res.txt
        timeout 1h ./main PRE $data 10000 2 $p 6 160 >> res.txt
    done

    for l in 2 4 6 8 10
    do
        timeout 1h ./main BSL $data 10000 2 0.5 $l >> res.txt
        timeout 1h ./main DP $data 10000 2 0.5 $l >> res.txt
        timeout 1h ./main DP+ $data 10000 2 0.5 $l >> res.txt
        timeout 1h ./main GDY $data 10000 2 0.5 $l 0.01 >> res.txt
        timeout 1h ./main PDBB $data 10000 2 0.5 $l >> res.txt
        timeout 1h ./main PRE $data 10000 2 0.5 $l 160 >> res.txt
    done
done
