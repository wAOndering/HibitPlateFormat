import glob
import os
import pandas as pd
import glob

print('')
print('-------------------------------------------')
print('How many plate needs to be averageds:')

bal = input('drag the file here')
bal = bal.replace('\\','/')
bal = bal.replace('"','')
print(bal)
print(glob.glob(bal+os.sep+'*.csv'))
# a = pd.read_csv(bal)

