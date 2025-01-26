import csv
import pandas as pd
import random
import numpy as np

m = 0
T = []
M = []
D = []
with open('iip_98_2000.txt', 'r') as fp, open('iip.dat', 'w', encoding='UTF-8') as fp_data:
    lines = fp.readlines()
    # print(lines[0].split())
    for line in lines:
        line = line.split()
        if (len(line) != 20): continue
        if (line[6] != 'R/V' and line[6] != 'VIS' and line[6] != 'RAD') : continue
        m += 1
        type = line[6]
        T.append(type)
        melt = float(line[16])
        days = float(line[17])
        M.append(melt + random.uniform(0, 1))
        D.append(days + random.uniform(0, 1))
    print(m)

    for i in range(len(T)):
        print(i, M[i], D[i], T[i])
        fp_data.write(str(M[i]) + ' ')
        fp_data.write(str(D[i]) + ' ')
        if T[i] == 'R/V':
            fp_data.write(str(0.8) + '\n')
        elif T[i] == 'VIS':
            fp_data.write(str(0.7) + '\n')
        elif T[i] == 'RAD':
            fp_data.write(str(0.6) + '\n')
