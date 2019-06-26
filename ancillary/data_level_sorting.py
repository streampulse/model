import sys
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from datetime import datetime, timedelta
import pandas as pd
import os
import re

app_dir = '/home/aaron/sp'
# app_dir = '/home/mike/git/streampulse/server_copy/sp'
sys.path.insert(0, app_dir)
os.chdir(app_dir)
import config as cfg

app = Flask(__name__)

app.config['SECRET_KEY'] = cfg.SECRET_KEY
app.config['SQLALCHEMY_DATABASE_URI'] = cfg.SQLALCHEMY_DATABASE_URI
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = cfg.SQLALCHEMY_TRACK_MODIFICATIONS
app.config['UPLOAD_FOLDER'] = cfg.UPLOAD_FOLDER
app.config['META_FOLDER'] = cfg.META_FOLDER
app.config['REACH_CHAR_FOLDER'] = cfg.REACH_CHAR_FOLDER
app.config['RESULTS_FOLDER'] = cfg.RESULTS_FOLDER
app.config['BULK_DNLD_FOLDER'] = cfg.BULK_DNLD_FOLDER
app.config['SECURITY_PASSWORD_SALT'] = cfg.SECURITY_PASSWORD_SALT

db = SQLAlchemy(app)


cols = pd.read_sql("select distinct region, site, concat(region, '_', site) as " +\
    "regsite, rawcol, dbcol from cols", db.engine)
cols['upload_file'] = ''
cols = cols.sort_values(['region', 'site'])

upfiles = os.listdir(app.config['UPLOAD_FOLDER'])
# i=0; j=0
for i in xrange(cols.shape[0]):
    sitefiles = [x for x in upfiles if cols.regsite[i] in x]
    for j in xrange(len(sitefiles)):
        header = os.popen('head ' + app.config['UPLOAD_FOLDER'] + '/' + sitefiles[j]).read()
        if cols.rawcol[i].encode('utf-8') in header:
            cols.upload_file[i] = sitefiles[j]
            continue

cols.to_csv('upfiles_by_col.csv')
