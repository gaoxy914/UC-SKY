#! /bin/bash

for data in ind/ anti/ corr/
do
    ./main PDBB $data 10000 6 0.5 6 0.01 >> res.txt
    for theta in 10 20 40 80 160 320 640 1280
    do
        ./main PRE $data 10000 6 0.5 6 $theta >> res.txt
    done
done

for data in ind/ anti/ corr/
do
    for vare in 0.1 0.01 0.001 0.0001 0.00001
    do
        ./main GDY $data 10000 6 0.5 6 $vae >> res.txt
    done
done

echo "n = 100000" >> res.txt

for data in ind/ anti/ corr/
do
    ./main PDBB $data 100000 6 0.5 6 0.01 >> res.txt
    for theta in 10 20 40 80 160 320 640 1280
    do
        ./main PRE $data 100000 6 0.5 6 $theta >> res.txt
    done
done

for data in ind/ anti/ corr/
do
    for vare in 0.1 0.01 0.001 0.0001 0.00001
    do
        ./main GDY $data 100000 6 0.5 6 $vae >> res.txt
    done
done

echo "l = 10" >> res.txt

for data in ind/ anti/ corr/
do
    ./main PDBB $data 100000 6 0.5 10 0.01 >> res.txt
    for theta in 10 20 40 80 160 320 640 1280
    do
        ./main PRE $data 100000 6 0.5 10 $theta >> res.txt
    done
done

for data in ind/ anti/ corr/
do
    for vare in 0.1 0.01 0.001 0.0001 0.00001
    do
        ./main GDY $data 100000 6 0.5 10 $vae >> res.txt
    done
done