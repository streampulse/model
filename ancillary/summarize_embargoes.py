import sys
sys.path.insert(0, '/home/mike/git/streampulse/server_copy/sp')

from flask import (Flask, Markup, session, flash, render_template, request,
    jsonify, url_for, make_response, send_file, redirect, g, send_from_directory)
from flask_login import LoginManager, login_user, logout_user, current_user, login_required
from werkzeug.security import generate_password_hash, check_password_hash
from sunrise_sunset import SunriseSunset as suns
from werkzeug.utils import secure_filename
from flask_sqlalchemy import SQLAlchemy
import sqlalchemy
from itsdangerous import URLSafeTimedSerializer
from datetime import datetime, timedelta
from dateutil import parser as dtparse
from math import log, sqrt, floor
import simplejson as json
from sklearn import svm
from operator import itemgetter
import pandas as pd
import numpy as np
import requests
import binascii
import tempfile
import zipfile
import shutil
import re
import os
import config as cfg
import readline #needed for rpy2 import in conda env
os.environ['R_HOME'] = '/usr/lib/R' #needed for rpy2 to find R. has to be a better way
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import markdown
import time
import traceback
import regex

app = Flask(__name__)

app.config['SECRET_KEY'] = cfg.SECRET_KEY
app.config['SQLALCHEMY_DATABASE_URI'] = cfg.SQLALCHEMY_DATABASE_URI
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = cfg.SQLALCHEMY_TRACK_MODIFICATIONS
app.config['UPLOAD_FOLDER'] = cfg.UPLOAD_FOLDER
app.config['META_FOLDER'] = cfg.META_FOLDER
app.config['GRAB_FOLDER'] = cfg.GRAB_FOLDER
app.config['RESULTS_FOLDER'] = cfg.RESULTS_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 700 * 1024 * 1024 # originally set to 16 MB; now 700
app.config['SECURITY_PASSWORD_SALT'] = cfg.SECURITY_PASSWORD_SALT

db = SQLAlchemy(app)
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = 'login'

#classes for SQLAlchemy's ORM
class Data(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    region = db.Column(db.String(10))
    site = db.Column(db.String(50))
    DateTime_UTC = db.Column(db.DateTime)
    variable = db.Column(db.String(50))
    value = db.Column(db.Float)
    flag = db.Column(db.Integer)
    upload_id = db.Column(db.Integer)

    def __init__(self, region, site, DateTime_UTC, variable, value, flag, upid):
        self.region = region
        self.site = site
        self.DateTime_UTC = DateTime_UTC
        self.variable = variable
        self.value = value
        self.flag = flag
        self.upload_id = upid

    def __repr__(self):
        return '<Data %r, %r, %r, %r, %r>' % (self.region, self.site,
        self.DateTime_UTC, self.variable, self.upload_id)

class Grabdata(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    region = db.Column(db.String(10))
    site = db.Column(db.String(50))
    DateTime_UTC = db.Column(db.DateTime)
    variable = db.Column(db.String(50))
    value = db.Column(db.Float)
    method = db.Column(db.String(40))
    write_in = db.Column(db.String(40))
    addtl = db.Column(db.String(40))
    flag = db.Column(db.Integer)
    upload_id = db.Column(db.Integer)

    def __init__(self, region, site, DateTime_UTC, variable, value, method,
        write_in, addtl, flag, upid):
        self.region = region
        self.site = site
        self.DateTime_UTC = DateTime_UTC
        self.variable = variable
        self.value = value
        self.method = method
        self.write_in = write_in
        self.addtl = addtl
        self.flag = flag
        self.upload_id = upid

    def __repr__(self):
        return '<Grabdata %r, %r, %r, %r, %r>' % (self.region, self.site,
            self.DateTime_UTC, self.variable, self.method, self.write_in,
            self.addtl, self.upload_id)

class Flag(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    region = db.Column(db.String(10))
    site = db.Column(db.String(50))
    startDate = db.Column(db.DateTime)
    endDate = db.Column(db.DateTime)
    variable = db.Column(db.String(50))
    flag = db.Column(db.String(50))
    comment = db.Column(db.String(255))
    by = db.Column(db.Integer) # user ID
    def __init__(self, region, site, startDate, endDate, variable, flag, comment, by):
        self.region = region
        self.site = site
        self.startDate = startDate
        self.endDate = endDate
        self.variable = variable
        self.flag = flag
        self.comment = comment
        self.by = by
    def __repr__(self):
        return '<Flag %r, %r, %r>' % (self.flag, self.comment, self.startDate)

class Grabflag(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    region = db.Column(db.String(10))
    site = db.Column(db.String(50))
    startDate = db.Column(db.DateTime)
    endDate = db.Column(db.DateTime)
    variable = db.Column(db.String(50))
    flag = db.Column(db.String(50))
    comment = db.Column(db.String(255))
    by = db.Column(db.Integer) # user ID
    def __init__(self, region, site, startDate, endDate, variable, flag, comment, by):
        self.region = region
        self.site = site
        self.startDate = startDate
        self.endDate = endDate
        self.variable = variable
        self.flag = flag
        self.comment = comment
        self.by = by
    def __repr__(self):
        return '<Grabflag %r, %r, %r>' % (self.flag, self.comment, self.startDate)

class Site(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    region = db.Column(db.String(10))
    site = db.Column(db.String(50))
    name = db.Column(db.String(100))
    latitude = db.Column(db.Float)
    longitude = db.Column(db.Float)
    usgs = db.Column(db.String(20))
    addDate = db.Column(db.DateTime)
    embargo = db.Column(db.Integer)
    by = db.Column(db.Integer)
    contact = db.Column(db.String(50))
    contactEmail = db.Column(db.String(255))
    def __init__(self, region, site, name, latitude, longitude, usgs, addDate, embargo, by, contact, contactEmail):
        self.region = region
        self.site = site
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.usgs = usgs
        self.addDate = addDate
        self.embargo = embargo
        self.by = by
        self.contact = contact
        self.contactEmail = contactEmail
    def __repr__(self):
        return '<Site %r, %r>' % (self.region, self.site)

class Cols(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    region = db.Column(db.String(10))
    site = db.Column(db.String(50))
    rawcol = db.Column(db.String(100))
    dbcol = db.Column(db.String(100))
    def __init__(self, region, site, rawcol, dbcol):
        self.region = region
        self.site = site
        self.rawcol = rawcol
        self.dbcol = dbcol
    def __repr__(self):
        return '<Cols %r, %r, %r>' % (self.region, self.site, self.dbcol)

class Grabcols(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    region = db.Column(db.String(10))
    site = db.Column(db.String(50))
    rawcol = db.Column(db.String(100))
    dbcol = db.Column(db.String(100))
    method = db.Column(db.String(40))
    write_in = db.Column(db.String(40))
    addtl = db.Column(db.String(40))

    def __init__(self, region, site, rawcol, dbcol, method, write_in, addtl):
        self.region = region
        self.site = site
        self.rawcol = rawcol
        self.dbcol = dbcol
        self.method = method
        self.write_in = write_in
        self.addtl = addtl

    def __repr__(self):
        return '<Grabcols %r, %r, %r>' % (self.region, self.site, self.dbcol)

class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(55), unique=True, index=True)
    password = db.Column(db.String(255))
    token = db.Column(db.String(100), nullable=False, server_default='')
    email = db.Column(db.String(255), unique=True)
    registered_on = db.Column(db.DateTime())
    confirmed = db.Column(db.Boolean)
    qaqc = db.Column(db.Text) # which qaqc sites can they save, comma separated?
    def __init__(self, username, password, email):
        self.username = username
        self.set_password(password)
        self.token = binascii.hexlify(os.urandom(10))
        self.email = email
        self.registered_on = datetime.utcnow()
        self.confirmed = True # do they agree to the policy?
        self.qaqc = ""
    def set_password(self, password):
        self.password = generate_password_hash(password)
    def check_password(self, password):
        return check_password_hash(self.password, password)
    def is_authenticated(self):
        return True
    def is_active(self):
        return True
    def is_anonymous(self):
        return False
    def get_id(self):
        return unicode(self.id)
    def qaqc_auth(self):
        return self.qaqc.split(",") # which tables can they edit
    def __repr__(self):
        return '<User %r>' % self.username

class Downloads(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    timestamp = db.Column(db.DateTime)
    userID = db.Column(db.Integer)
    email = db.Column(db.String(255))
    dnld_sites = db.Column(db.String(500))
    dnld_date0 = db.Column(db.DateTime)
    dnld_date1 = db.Column(db.DateTime)
    dnld_vars = db.Column(db.String(500))
    def __init__(self, timestamp, userID, email, dnld_sites, dnld_date0, dnld_date1, dnld_vars):
        self.timestamp = timestamp
        self.userID = userID
        self.email = email
        self.dnld_sites = dnld_sites
        self.dnld_date0 = dnld_date0
        self.dnld_date1 = dnld_date1
        self.dnld_vars = dnld_vars
    def __repr__(self):
        return '<Download %r>' % (self.id)

class Upload(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    filename = db.Column(db.String(100))
    # version = db.Column(db.Integer)

    def __init__(self, filename):
        self.filename = filename

    def __repr__(self):
        return '<Upload %r>' % (self.filename)

class Model(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    region = db.Column(db.String(2))
    site = db.Column(db.String(50))
    start_date = db.Column(db.DateTime)
    end_date = db.Column(db.DateTime)
    requested_variables = db.Column(db.String(200))
    year = db.Column(db.Integer)
    run_finished = db.Column(db.DateTime)
    model = db.Column(db.String(17))
    method = db.Column(db.String(5))
    engine = db.Column(db.String(4))
    rm_flagged = db.Column(db.String(35))
    used_rating_curve = db.Column(db.Boolean)
    pool = db.Column(db.String(7))
    proc_err = db.Column(db.Boolean)
    obs_err = db.Column(db.Boolean)
    proc_acor = db.Column(db.Boolean)
    ode_method = db.Column(db.String(9))
    deficit_src = db.Column(db.String(13))
    interv = db.Column(db.String(12))
    fillgaps = db.Column(db.String(13))
    estimate_areal_depth = db.Column(db.Boolean)
    O2_GOF = db.Column(db.Float)
    GPP_95CI = db.Column(db.Float)
    ER_95CI = db.Column(db.Float)
    prop_pos_ER = db.Column(db.Float)
    prop_neg_GPP = db.Column(db.Float)
    ER_K600_cor = db.Column(db.Float)
    coverage = db.Column(db.Integer)
    kmax = db.Column(db.Float)
    current_best = db.Column(db.Boolean)

    def __init__(self, region, site, start_date, end_date,
        requested_variables, year, run_finished, model, method, engine,
        rm_flagged, used_rating_curve, pool, proc_err, obs_err, proc_acor,
        ode_method, deficit_src, interv, fillgaps, estimate_areal_depth,
        O2_GOF, GPP_95CI, ER_95CI, prop_pos_ER, prop_neg_GPP, ER_K600_cor,
        coverage, kmax, current_best):

        self.region = region
        self.site = site
        self.start_date = start_date
        self.end_date = end_date
        self.requested_variables = requested_variables
        self.year = year
        self.run_finished = run_finished
        self.model = model
        self.method = method
        self.engine = engine
        self.rm_flagged = rm_flagged
        self.used_rating_curve = used_rating_curve
        self.pool = pool
        self.proc_err = proc_err
        self.obs_err = obs_err
        self.proc_acor = proc_acor
        self.ode_method = ode_method
        self.deficit_src = deficit_src
        self.interv = interv
        self.fillgaps = fillgaps
        self.estimate_areal_depth = estimate_areal_depth
        self.O2_GOF = O2_GOF
        self.GPP_95CI = GPP_95CI
        self.ER_95CI = ER_95CI
        self.prop_pos_ER = prop_pos_ER
        self.prop_neg_GPP = prop_neg_GPP
        self.ER_K600_cor = ER_K600_cor
        self.coverage = coverage
        self.kmax = kmax
        self.current_best = current_best

    def __repr__(self):
        return '<Data %r, %r, %r, %r>' % (self.region, self.site, self.year, self.current_best)

class Grabupload(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    filename = db.Column(db.String(100))
    # version = db.Column(db.Integer)

    def __init__(self, filename):
        self.filename = filename

    def __repr__(self):
        return '<Grabupload %r>' % (self.filename)

class Grdo(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50))
    email = db.Column(db.String(50))
    addDate = db.Column(db.DateTime)
    embargo = db.Column(db.Integer)
    notes = db.Column(db.String(5000))
    dataFiles = db.Column(db.String(5000))
    metaFiles = db.Column(db.String(5000))

    def __init__(self, name, email, addDate, embargo, notes, dataFiles, metaFiles):
        self.name = name
        self.email = email
        self.addDate = addDate
        self.embargo = embargo
        self.notes = notes
        self.dataFiles = dataFiles
        self.metaFiles = metaFiles

    def __repr__(self):
        return '<Data %r, %r, %r>' % (self.name, self.email, self.addDate)

db.create_all()

# end mostly unnecessary boilerplate

sitedata = pd.read_sql('select region as Region, site as Site, name as ' +\
    'Name, latitude as Lat, longitude as Lon, contact as Contact, ' +\
    'contactEmail as Email, usgs as `USGS gage`, `by`,' +\
    'embargo as EmbargoDays, addDate as AddDate, variableList as Variables, ' +\
    'firstRecord, lastRecord from site;', db.engine)

sitecounts = pd.read_sql("select region as Region, site as Site," +\
    " count(*) as DataCount from data group by Region, Site", db.engine)

sitedata = pd.merge(sitedata, sitecounts, on=['Region', 'Site'], how='left')

#calculate remaining embargo days
timedeltas = datetime.utcnow() - sitedata.AddDate
days_past = timedeltas.map(lambda x: int(x.total_seconds() / 60 / 60 / 24))
sitedata['EmbargoDays'] = sitedata['EmbargoDays'] * 365 - days_past
sitedata.loc[sitedata['EmbargoDays'] <= 0, 'EmbargoDays'] = 0

#format varlist and date range
pd.set_option('display.max_colwidth', 500)
core_variables = ['DO_mgL', 'satDO_mgL', 'DOsat_pct', 'WaterTemp_C',
    'Depth_m', 'Level_m', 'Discharge_m3s', 'Light_PAR', 'Light_lux']

varcells = []
for x in sitedata.Variables:
    if x is None:
        varcells.append(x)
    else:
        var_arr = np.asarray(x.split(','))
        isCore = np.in1d(var_arr, core_variables)
        core = var_arr[isCore]
        not_core = var_arr[~isCore]
        if any(core):
            core = core[np.argsort(pd.match(core, core_variables))]
        not_core.sort()
        var_arr = ', '.join(np.concatenate((core, not_core)))
        varcells.append(var_arr)

for i in xrange(len(varcells)):
    if varcells[i] is None:
        varcells[i] = '-'

sitedata.Variables = varcells
fr = sitedata['firstRecord'].dt.strftime('%Y-%m-%d')
lr = sitedata['lastRecord'].dt.strftime('%Y-%m-%d')
timerange = fr + ' to ' + lr
sitedata['Coverage'] = timerange.apply(lambda x: x if x != 'NaT to NaT' else '-')

#sort, filter, export
sitedata = sitedata.fillna('-').sort_values(['EmbargoDays', 'Region', 'Site'],
    ascending=False)
sitedata = sitedata.loc[~sitedata.by.isin([-900, -902]),]
sitedata = sitedata.drop(['Lat', 'Lon', 'USGS gage', 'by', 'firstRecord',
    'lastRecord'], axis=1)
sitedata = sitedata[['Region','Site','Name','EmbargoDays','DataCount','Coverage',
    'AddDate','Contact','Email','Variables']]
sitedata.to_csv('~/Desktop/embargo_summary.csv', index=False, encoding='utf-8')
