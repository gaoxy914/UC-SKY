#! /bin/bash

for data in ind/ anti/ corr/
do
    for n in 100 1000 10000 50000 100000 500000 1000000
    do
        timeout 1h ./main BSL $data $n 6 0.5 6 >> res.txt
        timeout 1h ./main GDY $data $n 6 0.5 6 0.01 >> res.txt
        timeout 1h ./main PDBB $data $n 6 0.5 6 >> res.txt
        timeout 1h ./main PRE $data $n 6 0.5 6 160 >> res.txt
    done

    for d in 2 4 6 8 10
    do
        timeout 1h ./main BSL $data 10000 $d 0.5 6 >> res.txt
        timeout 1h ./main GDY $data 10000 $d 0.5 6 0.01 >> res.txt
        timeout 1h ./main PDBB $data 10000 $d 6 0.5 6 >> res.txt
        timeout 1h ./main PRE $data 10000 $d 0.5 6 160 >> res.txt
    done

    for p in 0.2 0.35 0.5 0.65 0.8
    do
        timeout 1h ./main BSL $data 10000 6 $p 6 >> res.txt
        timeout 1h ./main GDY $data 10000 6 $p 6 0.01 >> res.txt
        timeout 1h ./main PDBB $data 10000 6 $p 6 >> res.txt
        timeout 1h ./main PRE $data 10000 6 $p 6 160 >> res.txt
    done

    for l in 2 4 6 8 10
    do
        timeout 1h ./main BSL $data 10000 6 0.5 $l >> res.txt
        timeout 1h ./main GDY $data 10000 6 0.5 $l 0.01 >> res.txt
        timeout 1h ./main PDBB $data 10000 6 0.5 $l >> res.txt
        timeout 1h ./main PRE $data 10000 6 0.5 $l 160 >> res.txt
    done

    for e in 0 0.1 0.01 0.001 0.0001
    do
        timeout 1h ./main GDY $data 10000 6 0.5 6 $e >> res.txt
    done

    for theta in 10 20 40 80 160 320 640 1280
    do
        timeout 1h ./main PRE $data 10000 6 0.5 6 $theta >> res.txt
    done
done
