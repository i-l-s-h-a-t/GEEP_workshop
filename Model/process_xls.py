import pandas as pd

data = pd.read_excel('time_data.xlsx')

# list the column names that should be removed
press_gridcells = data.filter(like='reservoir').columns.tolist()
chem_cols = data.filter(like='Kmol').columns.tolist()

# remove columns from data
data.drop(columns=press_gridcells + chem_cols, inplace=True)

# add time in years
data['Time (years)'] = data['time'] / 365.25

x = data.columns.tolist()

for i in x:
    print(i)

