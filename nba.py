from email import header
from operator import index
# from nba_api.stats.static import players
# from nba_api.stats.endpoints import playercareerstats
# from nba_api.stats.endpoints import commonplayerinfo
import pandas as pd
import numpy as np
import os
import random


# nba_players = players.get_players()
# print('Number of players fetched: {}'.format(len(nba_players)))
# for player in nba_players:
    # career = playercareerstats.PlayerCareerStats(player_id=player['id'])
    # df = career.get_data_frames()[0]
    # if not os.path.exists('player.csv'):
        # df.to_csv('player.csv', mode='a', index=False, index_label=False)
    # else:
        # df.to_csv('player.csv', mode='a', index=False, index_label=False, header=False)


cols = ['MIN', 'FGM', 'PTS', 'TOV', 'BLK', 'STL', 'AST', 'REB']
df = pd.read_csv('player.csv')[cols]
df = df.dropna()
max_list = df.max().to_list()
print(max_list)
with open('nba.dat', 'w', encoding='UTF-8') as fp:
    print(len(df))
    # fp.write(str(len(df)) + '\n')
    for index, row in df.iterrows():
        fp.write(str(max_list[0] - row['MIN']) + ' ')
        fp.write(str(max_list[1] - row['FGM']) + ' ')
        fp.write(str(max_list[2] - row['PTS']) + ' ')
        fp.write(str(row['TOV']) + ' ')
        fp.write(str(max_list[4] - row['BLK']) + ' ')
        fp.write(str(max_list[5] - row['STL']) + ' ')
        fp.write(str(max_list[6] - row['AST']) + ' ')
        fp.write(str(max_list[7] - row['REB']) + ' ')
        prob = random.random()
        fp.write(str(round(prob, 6)) + '\n')