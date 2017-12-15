import os
import pandas as pd

#script for generating gradient boosted classification tree model to predict
#problematic datapoints in data uploaded to StreamPULSE data portal

#first downoad updated data and flags from https://data.streampulse.org/download

wd = '/home/mike/Dropbox/streampulse/data/NC_download/'
site_files = os.listdir(wd)

spsite = pd.read_csv(wd+site_files[0])
sitevars = spsite.variable.unique() #variables available for site
for i in sitevars:
spsite[spsite.variable == 'light_lux']
