import pandas as pd
import datetime

logpath = '/storage/home/len56/work/warm_jupiters/data/datalog.csv'
df = pd.DataFrame(columns=['tag'])
df.to_csv(logpath,index=False)