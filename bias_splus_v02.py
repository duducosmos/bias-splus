#! /usr/bin python
# -*- coding: utf-8 -*-

"""
Name: bias_splus
Version 0.3
Description: Provide statistical analysis for s-plus bias images over periods of time
Author: Walter Santos
Created on: Set 12th 2016
Last Updated: Nov 26th 2016
Latest Changes:
- Include the acceptable percentage (above which the bias is flagged
as not acceptable) into a command line parameter
- Changed default output filename to 'bias_splus_$finalDate'
- Include errors and warning in he output file log
- Change the way we look for new BIAS images using a mongo db setupt by Tiago/William
Instructions: Run the code with a -h/--help option to list all the arguments,
necessary and optional, on the command line.
Requirements: Numpy, Astropy, Matplotlib
"""

"""
NOTE:
For human input, unless stated otherwise, date should be given as a string by yyyy/mm/dd, 
for example, '2016/09/09'
For function input, are given as python date objects
"""

DEFAULT_INITIAL_DATE = '2016/09/01'
#DEFAULT_FINAL_DATE = '2016/09/11'
DEFAULT_ACCEPTABLE_PERC = 0.05


"""
IMPORTS
"""
import matplotlib as mpl
mpl.use('Agg')
import argparse
import os
import os.path
import numpy as np
from astropy.io import fits
from datetime import date, timedelta
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib.dates import DayLocator, DateFormatter, date2num
from pymongo import MongoClient
import time

"""
FUNCTIONS
"""
def outputPlot(fulldate, xdates, ypoints, errors, outfile, perc):
    fig, ax = plt.subplots()

    y = np.asarray(ypoints,dtype=np.float)
    err = np.asarray(errors)
    x = np.asarray(date2num(xdates))

    colors = ['red', 'blue']
    levels = [0, 1]

    cmap, norm = mpl.colors.from_levels_and_colors(levels=levels, colors=colors, extend='max')

    colorredalert = np.where((y >= (1.0+perc)) | (y <= (1.0-perc)), 0, 1)

    ax.scatter(x,y,c=colorredalert, marker='o', edgecolor='none', cmap=cmap, norm=norm, zorder=1)   

    ax.errorbar(x, y, yerr=err, fmt='none', ecolor='k', zorder=0) 
 
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.set_ylabel('Counts normalized by master') #Ylabel
    ax.xaxis.set_major_locator(DayLocator())
    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
    ax.fmt_xdata = DateFormatter('%Y-%m-%d')
    ax.set_title('Mean Bias Variation')
    ax.set_xlim(fulldate[0]-timedelta(days=1), fulldate[-1]+timedelta(days=1))
    ax.set_ylim(min(y)-5*max(err), max(y)+5*max(err))
    fig.autofmt_xdate()    
    
    plt.savefig(outfile,dpi=80) # Saving figure


def daterange(start_date, end_date):
    for ordinal in range(start_date.toordinal(), end_date.toordinal()):
        yield date.fromordinal(ordinal)

def getDatesLists(init, final):
    list_dates = daterange(init, final)
    return [d for d in list_dates], [str(d.year)+('%02d' % d.month)+('%02d' % d.day) for d  in list_dates]



"""
MAIN
"""
parser = argparse.ArgumentParser(description='S-PLUS Bias Analysis')
parser.add_argument('--initdate', default=None, help='initial \
                    date for the analysis', metavar='yyyy/mm/dd')
parser.add_argument('--finaldate', default=None, help='final \
                    date for the analysis', metavar='yyyy/mm/dd')
parser.add_argument('-o','--output', default=None, help='output file name, also \
                     used for the output plot name, an appropriate default \
                     will be given if None', metavar='output_file')
parser.add_argument('-m','--master', default='splus_master_bias.fits', help='input \
                     master bias image', metavar='MASTER_BIAS.fits')
parser.add_argument('--perc', default=DEFAULT_ACCEPTABLE_PERC, type=float,
                    help='percentage above which the bias is flagged \
                    as not acceptable')
parser.add_argument('--host', help='The host of the mongo database. \
                    Mandatory if not ', type=str, default='192.168.20.118')
parser.add_argument('--port', help='The port of the mongo database. \
                    Mandatory if not ', type=str, default='27017')
parser.add_argument('--path', help='path to the images', type=str, default='27017')

#command line example:
"""
command_line = '--initdate "2016/09/05" \
                --finaldate "2016/09/11" \
                -o "bias_splus_05_10" \
                -m "splus_master_bias.fits"

python bias_splus_v02.py --initdate "2016/09/05" --finaldate "2016/09/11" -o "bias_splus_05_10" -m "splus_master_bias.fits"
or for new mode, ony takes into account new bias/dates:
python bias_splus_v02.py -o "bias_splus_26_sept" -m "splus_master_bias.fits"
"""

"""
Parse the input parameters, define and declare variables
"""
args = parser.parse_args()
settings = vars(args)

host = settings['host']
port = settings['port']
path = settings['path']

client = MongoClient('%s:%s' % (host, port))

if ((settings['initdate'] and settings['finaldate']) is not None):
    initDate = date(int(settings['initdate'].split('/')[0]),int(settings['initdate'].split('/')[1]),int(settings['initdate'].split('/')[2]))
    finalDate = date(int(settings['finaldate'].split('/')[0]),int(settings['finaldate'].split('/')[1]),int(settings['finaldate'].split('/')[2]))
    fulldatelist, datelist = getDatesLists(initDate,finalDate)
else:
    initDate = date(int(DEFAULT_INITIAL_DATE.split('/')[0]),int(DEFAULT_INITIAL_DATE.split('/')[1]),int(DEFAULT_INITIAL_DATE.split('/')[2]))
    finalDate = date.today()
    fulldatelist, datelist = getDatesLists(initDate,finalDate)
    
if(settings['output']==None):
    outputFilename = 'bias_splus'+str(finalDate)
else:
    outputFilename = settings['output']
fout = open(outputFilename+'.log', "w")

masterFilename = settings['master']
try:
    _mdata = fits.open(masterFilename)[0].data
    m_nx, m_ny = np.shape(_mdata)
    mdata = np.zeros((m_nx, m_ny),dtype=np.float32)
    mdata += _mdata
except:
    print('ERROR: Could not find or access data of master bias named: '+masterFilename)
    fout.write('\n ERROR: Could not find or access data of master bias named: '+masterFilename+'\n')
    raise

acceptablePerc = settings['perc']

bias_day_mean = list()
bias_day_error = list()
actualdates = datelist

for night in datelist:
    cursor = client.images.fits_keywords.find({'_night': night,'IMAGETYP': 'ZERO'}).sort([('OBJECT', 1), ('FILENAME', 1), ])
    nbias = cursor.count()
    if nbias <= 0:
        print('WARNING: There are no suitable bias files for date '+night+'!')
        fout.write('\n WARNING: There are no suitable bias files for date '+night+'!')
        actualdates.remove(night)
    else:
        fout.write('\n\n# All bias for date '+ night+' ...')
        fout.write('\n# date filename mean median max min std mean_per_master median_per_master std_per_master')
        print '\n\n# All bias for date '+ night+' ...'
        print '# date filename mean median max min std mean_per_master median_per_master std_per_master'
    
        biasfilelist = [str(frame['FILENAME']) for frame in cursor]
        bias_stack = np.zeros((nbias,m_nx,m_ny),dtype=np.float32)
        bias_stack_per_master = np.zeros((nbias,m_nx,m_ny),dtype=np.float32)
        i = 0
        for b in biasfilelist:
            try:
                _bdata = fits.open(path+'/'+night+'/'+b)[0].data
                nx, ny = np.shape(_bdata)
                if (nx != m_nx or ny != m_ny):
                    raise		
                bdata = np.zeros((nx,ny),dtype=np.float32)
                bdata += _bdata
            except:
                print('WARNING: Unable to open or corrupted file '+b+'. Ignoring it.\n')
                fout.write('\n WARNING: Unable to open or corrupted file '+b+'. Ignoring it.\n')
                bias_stack = np.delete(bias_stack,-1,axis=0)
                nbias-=1
            else:
                data_per_master = bdata/mdata
                print night, b, int(np.mean(bdata)), int(np.median(bdata)), int(np.max(bdata)), int(np.min(bdata)), np.std(bdata), np.mean(data_per_master), np.median(data_per_master), np.std(data_per_master)
                fout.write('\n'+night+' '+b+' '+str(int(np.mean(bdata)))+' '+str(int(np.median(bdata)))+' '+str(int(np.max(bdata)))+' '+str(int(np.min(bdata)))+' '+str(np.std(bdata))+' '+str(np.mean(data_per_master))+' '+str(np.median(data_per_master))+' '+str(np.std(data_per_master)))
                bias_stack[i] = bdata
                bias_stack_per_master[i] = data_per_master
                i+=1
    fout.write('\n# Combined ' + str(nbias) + 'bias for the date '+ night+' ...')
    fout.write('\n# date mean median max min std mean_per_master median_per_master std_per_master')
    print '# Combined ' + str(nbias) + ' bias for the date '+ night+' ...'
    print '# date mean median max min std mean_per_master median_per_master std_per_master'
    stacked_bias = np.median(bias_stack,axis=0)
    stacked_bias_per_master = stacked_bias/mdata
    fout.write('\n'+night+' '+str(int(np.mean(stacked_bias)))+' '+str(int(np.median(stacked_bias)))+' '+ str(int(np.max(stacked_bias)))+' '+ str(int(np.min(stacked_bias)))+' '+ str(np.std(stacked_bias))+' '+str(np.mean(stacked_bias_per_master))+' '+str(np.median(stacked_bias_per_master))+' '+str(np.std(stacked_bias_per_master)))
    print night, int(np.mean(stacked_bias)), int(np.median(stacked_bias)), int(np.max(stacked_bias)), int(np.min(stacked_bias)), np.std(stacked_bias), np.mean(stacked_bias_per_master), np.median(stacked_bias_per_master), np.std(stacked_bias_per_master)
    if((np.mean(stacked_bias_per_master) >= (1.0+acceptablePerc)) or (np.mean(stacked_bias_per_master) <= (1.0-acceptablePerc))):
        print('WARNING: Bias for the date '+ night+'above acceptable value compared to master!')
        fout.write('\n WARNING: Bias for the date '+ night+'above acceptable value compared to master! \n')
    bias_day_mean.append(np.mean(stacked_bias_per_master))
    bias_day_error.append(np.std(stacked_bias_per_master))
    

outputPlot(fulldatelist, actualdates, bias_day_mean, bias_day_error, settings['output']+'.png', acceptablePerc)
fout.close()