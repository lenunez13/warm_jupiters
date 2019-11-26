import pandas as pd
import datetime

date = str(datetime.datetime.now()).split()[0]
df = pd.DataFrame({'date':[date]})
df.to_csv('/storage/home/len56/work/warm_jupiters/data/datalog.csv',index=False)