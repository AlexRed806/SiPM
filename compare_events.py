import pandas as pd

table_1 = pd.read_csv("/home/admin/Desktop/SiPM/results/HPK_coincidences/uncovered/horizontal_4ov_c1.csv")
table_2 = pd.read_csv("/home/admin/Desktop/SiPM/results/HPK_coincidences/uncovered/horizontal_4ov_c2.csv")

n_coincidences = 0
coinc = 0

for x in table_1['Timestamp']:
        coinc = len(table_2[ abs(table_2['Timestamp'] - x) < 0.0000001])
        #print(coinc)
        if coinc == 1: n_coincidences += 1
        if coinc >1: print ("Too many")

print(n_coincidences)
