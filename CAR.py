import pandas as pd
import random

cols = ['price', 'yearOfRegistration', 'powerPS', 'kilometer', 'dateCreated']
df = pd.read_csv('autos.csv')[cols]
df = df.dropna()
df = df[df['price'] != 0]
df = df[df['price'] < 100000]
df = df[df['powerPS'] != 0]
df = df[df['yearOfRegistration'] < 2023]
max_list = df.max().to_list()
min_list = df.min().to_list()
print(max_list)
print(min_list)
with open('car.dat', 'w', encoding='UTF-8') as fp:
    print(len(df))
#    fp.write(str(len(df)) + ' ')
    for index, row in df.iterrows():
        fp.write(str(row['price']) + ' ')
        fp.write(str(max_list[1] - row['yearOfRegistration']) + ' ')
        fp.write(str(max_list[2] - row['powerPS']) + ' ')
        fp.write(str(row['kilometer']) + ' ')
        
        date = row['dateCreated']
        date = date.split()[0]
        day = int(date.split('/')[2])
        car_prob = 1
        sold_prob = random.random()
        sold = False
        while not sold and day < 42:
            if random.random() > sold_prob:
                car_prob *= (1 - sold_prob)
            else:
                car_prob *= sold_prob
                sold = True
            day += 1
        fp.write(str(round(car_prob, 6)) + '\n')