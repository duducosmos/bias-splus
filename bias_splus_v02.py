#! /usr/bin python
# -*- coding: utf-8 -*-

"""
Name: bias_splus
Version 0.2
Description: Provide statistical analysis for s-plus bias images over periods of time
Author: Walter Santos
Created on: Set 12th 2016
Last Updated: Nov 26th 2016
Latest Changes:
- Try/exception clauses to ignore bad bias files
- Added a new input parameter: master bias
- Results include biases division by this master, including thefinal plot
- When no dates are given, only new downaloaded (updated in the server
    are considered for the analysis)
- If a bias for a given night is off by a certain percentage compared to the master,
    a warning is printed out and the point gets redin the plot, instead of the default blue
Instructions: Run the code with a -h/--help option to list all the arguments,
necessary and optional, on the command line.
Requirements: Numpy, Astropy, Matplotlib

----
EXAMPLE:

"""

"""
NOTE:
For human input, unless stated otherwise, date should be given as a string by yyyy/mm/dd, 
for example, '2016/09/09'
For function input, are given as python date objects
"""

DEFAULT_INITIAL_DATE = '2016/09/01'
#DEFAULT_FINAL_DATE = '2016/09/11'
ACCEPTABLE_PERC = 0.05


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

"""
FUNCTIONS
"""
def outputPlot(fulldate, xdates, ypoints, errors, outfile):
    fig, ax = plt.subplots()

    y = np.asarray(ypoints,dtype=np.float)
    err = np.asarray(errors)
    x = np.asarray(date2num(xdates))

    colors = ['red', 'blue']
    levels = [0, 1]

    cmap, norm = mpl.colors.from_levels_and_colors(levels=levels, colors=colors, extend='max')

    colorredalert = np.where((y >= (1.0+ACCEPTABLE_PERC)) | (y <= (1.0-ACCEPTABLE_PERC)), 0, 1)

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
        

def getFilelistsDates(init, final):
    list_dates = daterange(init, final)
    filelist = ['filelist_'+str(d.year)+('%02d' % d.month)+('%02d' % d.day)+'.csv' for d in list_dates]
    return filelist
    
    
def downloadFilelistsCSV(filelist):
    csvlist = list()
    cmd = 'wget --no-check-certificate -c https://t80s_images:t80s_images_keywords_pass@splus.astro.ufsc.br/'

    for f in filelist:
        if os.path.isfile(f):
            csvlist.append(f)
        else:
            cmd_wget = cmd+f
            os.system(cmd_wget)
            if os.path.isfile(f):
                csvlist.append(f)
    return csvlist


def downloadNewFilelistsCSV(filelist):
    csvlist = list()
    cmd = 'wget --no-check-certificate -c https://t80s_images:t80s_images_keywords_pass@splus.astro.ufsc.br/'
    
    for f in filelist:
        if os.path.isfile(f) == False:
            cmd_wget = cmd+f
            os.system(cmd_wget)
            if os.path.isfile(f):
                csvlist.append(f)
                
    return csvlist        
        
    
def getBiaslistFromCSV(csvfile):
    biaslist = list()
    biasnames = list()
    f = open(csvfile, 'r')
    for line in f:
        line = line.strip()
        columns = line.split()
        pathfile = columns[2]
        namefile = columns[2].split('/')[-1]
        for c in columns[3].split(','):
            if c == 'ZERO':
                biaslist.append(pathfile)
                biasnames.append(namefile)
    return biaslist, biasnames
        
        
def convertCSVtoDate(csvfile):
    d = csvfile[9:-4]
    return date(int(d[0:4]), int(d[4:6]), int(d[6:8]))

    
def downloadBiasFromDate(biaslist):
    cmd = 'wget --no-check-certificate -c '
    for b in biaslist:
        namefile = b.split('/')[-1]
        if os.path.isfile(namefile) == False:
            cmd_get = cmd+b
            os.system(cmd_get)
    
"""
MAIN
"""
parser = argparse.ArgumentParser(description='S-PLUS Bias Analysis')
parser.add_argument('--initdate', default=None, help='initial \
                    date for the analysis', metavar='yyyy/mm/dd')
parser.add_argument('--finaldate', default=None, help='final \
                    date for the analysis', metavar='yyyy/mm/dd')
parser.add_argument('-o','--output', default='bias_splus_sept', help='output file name, also \
                     used for the output plot name, an appropriate default \
                     will be given if None', metavar='output_file')
parser.add_argument('-m','--master', default='splus_master_bias.fits', help='input \
                     master bias image', metavar='MASTER_BIAS.fits')

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

if ((settings['initdate'] and settings['finaldate']) != None):
    initDate = date(int(settings['initdate'].split('/')[0]),int(settings['initdate'].split('/')[1]),int(settings['initdate'].split('/')[2]))
    finalDate = date(int(settings['finaldate'].split('/')[0]),int(settings['finaldate'].split('/')[1]),int(settings['finaldate'].split('/')[2]))
    list_dates = daterange(initDate, finalDate)
    datelist = getFilelistsDates(initDate,finalDate)
    csvlist = downloadFilelistsCSV(datelist)
    fulldatelist = [d for d in list_dates]
    actualdates = map(convertCSVtoDate, csvlist)
else:
    initDate = date(int(DEFAULT_INITIAL_DATE.split('/')[0]),int(DEFAULT_INITIAL_DATE.split('/')[1]),int(DEFAULT_INITIAL_DATE.split('/')[2]))
    finalDate = date.today()
    datelist = getFilelistsDates(initDate,finalDate)
    csvlist = downloadNewFilelistsCSV(datelist)
    actualdates = map(convertCSVtoDate, csvlist)
    fulldatelist = actualdates

masterFilename = settings['master']
try:
    _mdata = fits.open(masterFilename)[0].data
    m_nx, m_ny = np.shape(_mdata)
    mdata = np.zeros((m_nx, m_ny),dtype=np.float32)
    mdata += _mdata
except:
    print('Could not find or access data of master bias named: '+masterFilename)
    raise
    
fout = open(settings['output']+'.dat', "w")

bias_day_mean = list()
bias_day_error = list()

for csv in csvlist:
    dateback = convertCSVtoDate(csv)
    fout.write('\n\n# All bias for date '+ str(dateback)+' ...')
    fout.write('\n# date filename mean median max min std mean_per_master median_per_master std_per_master')
    print '\n\n# All bias for date '+ str(dateback)+' ...'
    print '# date filename mean median max min std mean_per_master median_per_master std_per_master'
    biasfilelist, biasnameslist = getBiaslistFromCSV(csv)
    downloadBiasFromDate(biasfilelist)
    
    nbias = len(biasnameslist)
   
    bias_stack = np.zeros((nbias,m_nx,m_ny),dtype=np.float32)
    bias_stack_per_master = np.zeros((nbias,m_nx,m_ny),dtype=np.float32)
    
    i = 0
    for b in biasnameslist:
        try:
            _bdata = fits.open(b)[0].data
            nx, ny = np.shape(_bdata)
            if (nx != m_nx or ny != m_ny):
		raise		
            bdata = np.zeros((nx,ny),dtype=np.float32)
            bdata += _bdata
        except:
            print('WARNING: Unable to open or corrupted file '+b+'. Ignoring it.\n')
            bias_stack = np.delete(bias_stack,-1,axis=0)
            nbias-=1
        else:
            data_per_master = bdata/mdata
            print dateback, b, int(np.mean(bdata)), int(np.median(bdata)), int(np.max(bdata)), int(np.min(bdata)), np.std(bdata), np.mean(data_per_master), np.median(data_per_master), np.std(data_per_master)
            fout.write('\n'+str(dateback)+' '+b+' '+str(int(np.mean(bdata)))+' '+str(int(np.median(bdata)))+' '+str(int(np.max(bdata)))+' '+str(int(np.min(bdata)))+' '+str(np.std(bdata))+' '+str(np.mean(data_per_master))+' '+str(np.median(data_per_master))+' '+str(np.std(data_per_master)))
            bias_stack[i] = bdata
            bias_stack_per_master[i] = bdata/mdata
            i+=1
   
    if(nbias == 0):
        print('WARNING: There are no suitable bias files for date '+str(dateback)+'!')
        actualdates.remove(dateback)
    else:
        fout.write('\n# Combined ' + str(nbias) + 'bias for the date '+ str(dateback)+' ...')
        fout.write('\n# date mean median max min std mean_per_master median_per_master std_per_master')
        print '# Combined ' + str(nbias) + ' bias for the date '+ str(dateback)+' ...'
        print '# date mean median max min std mean_per_master median_per_master std_per_master'
        stacked_bias = np.median(bias_stack,axis=0)
        stacked_bias_per_master = stacked_bias/mdata
        fout.write('\n'+str(dateback)+' '+str(int(np.mean(stacked_bias)))+' '+str(int(np.median(stacked_bias)))+' '+ str(int(np.max(stacked_bias)))+' '+ str(int(np.min(stacked_bias)))+' '+ str(np.std(stacked_bias))+' '+str(np.mean(stacked_bias_per_master))+' '+str(np.median(stacked_bias_per_master))+' '+str(np.std(stacked_bias_per_master)))
        print dateback, int(np.mean(stacked_bias)), int(np.median(stacked_bias)), int(np.max(stacked_bias)), int(np.min(stacked_bias)), np.std(stacked_bias), np.mean(stacked_bias_per_master), np.median(stacked_bias_per_master), np.std(stacked_bias_per_master)
        if((np.mean(stacked_bias_per_master) >= (1.0+ACCEPTABLE_PERC)) or (np.mean(stacked_bias_per_master) <= (1.0-ACCEPTABLE_PERC))):
            print('WARNING: Bias for the date '+ str(dateback)+'above acceptable value compared to master!')
        bias_day_mean.append(np.mean(stacked_bias_per_master))
        bias_day_error.append(np.std(stacked_bias_per_master))
    
outputPlot(fulldatelist, actualdates, bias_day_mean, bias_day_error, settings['output']+'.png')
fout.close()