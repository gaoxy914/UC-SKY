# UC-SKY

## Compliation

Run makefile.

## Generate synthetic datasets

To generate all synthetic datasets, run ./main gen-data.

To generate 2D datasets with m groups, run ./main gen-data-2d m.

To generate nba.dat, run python3 nba.py.

To generate iip.dat, run python3 IIP.py.

To generate car.dat, run python3 CAR.py.

## Run algorithms

Commands format:
```
./main BSL/DP/DP+/GDY/SKY/PDBB data n d alpha l epsilon
```
|Parameters|Meaning|Example|
|:---:|---|:---:|
|data|data type|ind, anti, car, nba.data, car.dat, iip.dat|
|n|data size|10000, ...|
|d|data dimension|2, 3, 4, ...|
|alpha|probability center|0.2, ...|
|l|output size|1, 2, ...|
|epsilon|step for GDY|0.001, ...|

