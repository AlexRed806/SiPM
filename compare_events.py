import pandas as pd

table_1 = pd.read_csv()
table_2 = pd.read_cvs()

n_coincidences = 0

for x in table_1['Timestamp']:
	# these are all the similar ones
	if table_2[ abs(table_2['Timestamp'] - x) < 0.01]:
        n_coincidences += 1

print(n_coincidences)
