# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 07:45:50 2015

@author: Alex
"""

#timer
import datetime as dt
start_time = dt.datetime.now()
print 'Start time: '+start_time.strftime('%H:%M:%S')



#### Import modules
## Data Processing
import os
import numpy as np
import pandas as pd

## Statistical Analysis
from scipy import stats
import pandas.stats.moments as m
from scipy.stats import pearsonr as pearson_r
from scipy.stats import spearmanr as spearman_r

import statsmodels.formula.api as smf
import statsmodels.stats.api

## Call R modules
#from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com

## Make sure R is communicating
ro.r('x=c()')
ro.r('x[1]="lets talk to R"')
print(ro.r('x'))

import pypandoc

#### Plotting Tools
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec 
import pylab
from AnnoteFinder import AnnoteFinder
import seaborn as sns
plt.close('all')
plt.ion()

##custom modules
import misc_time
from misc_time import * 
from misc_numpy import *
from misc_matplotlib import * 

## Set Pandas display options
pd.set_option('display.large_repr', 'truncate')
pd.set_option('display.width', 180)
pd.set_option('display.max_rows', 30)
pd.set_option('display.max_columns', 13)

#### DIRECTORIES
git=True
if git==True: ## Git repository
    maindir = 'C:/Users/Alex/Documents/GitHub/Nuuuli-Sediment-Flux/' 
    datadir=maindir+'Data/'
    dataoutputdir = datadir+'Output/'
    GISdir = maindir+'Data/GIS/'
    figdir = maindir+'Figures/'
    tabledir = maindir+'Tables/'
    dirs={'main':maindir,'data':datadir,'GIS':GISdir,'fig':figdir}
elif git!=True: ## Local folders
    maindir = 'C:/Users/Alex/Desktop/'### samoa/
    csvoutputdir = datadir+'samoa/WATERSHED_ANALYSIS/FAGAALU/MasterDataFiles/csv_output/'
    savedir = datadir+'samoa/WATERSHED_ANALYSIS/GoodFigures/'
    figdir = datadir+'samoa/WATERSHED_ANALYSIS/GoodFigures/rawfigoutput/'

def show_plot(show=False,fig=figure):
    if show==True:
        plt.show()
def logaxes(log=False,fig=figure):
    if log==True:
        print 'log axes'
        for ax in fig.axes:
            ax.set_yscale('log'), ax.set_xscale('log')
    return
def savefig(save=True,filename=''):
    if save==True:
        #plt.savefig(filename+'.pdf') ## for publication
        plt.savefig(filename+'.png') ## for manuscript
    return
def pltdefault():
    global figdir
    plt.rcdefaults()
    #figdir = datadir+'samoa/WATERSHED_ANALYSIS/GoodFigures/rawfigoutput/'
    return  
    

def letter_subplots(fig,x=0.1,y=0.95,vertical='top',horizontal='right',Color='k',font_size=10,font_weight='bold'):
    sub_plot_count = 0
    sub_plot_letters = {0:'(a)',1:'(b)',2:'(c)',3:'(d)',4:'(e)',5:'(f)',6:'(g)',7:'(h)',8:'(i)'}
    for ax in fig.axes:
        ax.text(x,y,sub_plot_letters[sub_plot_count],verticalalignment=vertical, horizontalalignment=horizontal,transform=ax.transAxes,color=Color,fontsize=font_size,fontweight=font_weight)
        sub_plot_count+=1
    return 

    
## Figure formatting
#publishable =  plotsettings.Set('GlobEnvChange')    ## plotsettings.py
#publishable.set_figsize(n_columns = 2, n_rows = 2)
    
mpl.rc_file(maindir+'johrc.rc')
mpl.rcParams['savefig.directory']=maindir+'rawfig/'
mpl.rcParams
mpl.rc('legend',scatterpoints=1)  
## Ticks
my_locator = matplotlib.ticker.MaxNLocator(4)

    
def pltsns(style='white',context='paper'):
    global figdir
    sns.set_style(style)
    sns.set_style({'legend.frameon':True})
    sns.set_context(context)
    ## Some Formatting
    mpl.rcParams['savefig.directory']=maindir+'rawfig/'
    mpl.rcParams['savefig.format']= '.pdf'
    mpl.rcParams['figure.dpi']=100
    mpl.rcParams['savefig.dpi']= 300 # mpl.rcParams['figure.dpi']
    mpl.rcParams['axes.labelsize']= 11
    mpl.rcParams['xtick.labelsize']=11
    mpl.rcParams['ytick.labelsize']=11
    
    return
#pltsns()    

def xkcd():
    global figdir
    plt.xkcd()
    #figdir = datadir+'samoa/WATERSHED_ANALYSIS/GoodFigures/rawfigoutput/xkcd/'
    return

## Misc. plotting tools
def labelxy(i,x,y):
    annotes = pd.DataFrame([x,y]).apply(tuple,axis=0).values.tolist()
    annotestrings = ["%.1f"%an[0]+','+"%.1f"%an[1] for an in annotes]
    af = AnnoteFinder(x,y,annotestrings)
    pylab.connect('button_press_event', af)
    return

def labelindex(i,x,y,ax,display=False): 
    if ax==None:
        ax=plt.gca()
    indexstrings = [str(ind) for ind in i]
    if display ==True:
        for i in zip(indexstrings,x,y):
            print i
    af = AnnoteFinder(x,y,indexstrings,axis=ax)
    pylab.connect('button_press_event', af)
    return
    
def labelindex_subplot(ax,i,x,y): 
    indexstrings = [str(ind) for ind in i]
    af = AnnoteFinder(x,y,indexstrings,axis=ax)
    connect('button_press_event', af)
    return

def annotate_plot(frame,plot_col,label_col):
    frame = frame[frame[label_col].isnull()!=True]
    for label, x, y in zip(frame['fieldnotes'], frame.index, frame['SSC (mg/L)']):
            plt.annotate(label, xy=(x, y))
    return
    
def scaleSeries(series,new_scale=[100,10]):
    new_scale = new_scale
    OldRange = (series.max() - series.min())  
    NewRange = (new_scale[0] - new_scale[1])  
    NewSeriesValues = (((series - series.min()) * NewRange) / OldRange) + new_scale[1]
    return NewSeriesValues          
    
def power(x,a,b):
    y = a*(x**b)
    return y
    
def powerfunction(x,y,name='power rating',pvalue=0.01):
    ## put x and y in a dataframe so you can drop ones that don't match up  
    datadf = pd.DataFrame.from_dict({'x':x,'y':y}).dropna().apply(np.log10)   
    datadf = datadf[datadf>=-10] ##verify data is valid (not inf)
    regression = pd.ols(y=datadf['y'],x=datadf['x'])
    if pearson_r(datadf['x'],datadf['y'])[1] < pvalue:
        pearson = pearson_r(datadf['x'],datadf['y'])[0]
    else: 
        pearson = np.nan
    if  spearman_r(datadf['x'],datadf['y'])[1] < pvalue:
        spearman = spearman_r(datadf['x'],datadf['y'])[0]
    else:
        spearman = np.nan
    coeffdf = pd.DataFrame({'a':[10**regression.beta[1]],'b':[regression.beta[0]],
    'r2':[regression.r2],'rmse':[regression.rmse],'pearson':[pearson],'spearman':[spearman]},
index=[name])
    return coeffdf

def PowerFit(x,y,xspace = 'none',ax=plt,**kwargs):
    ## Develop power function for x and y
    powfunc = powerfunction(x,y) ## x and y should be Series
    a, b = powfunc['a'].values, powfunc['b'].values
    #print a,b
    if xspace == 'none':
        xvals = np.linspace(x.min()-10,x.max()*1.2)
        #print 'No xspace, calculating xvals: '+'%.0f'%x.min()+'-'+'%.0f'%x.max()+'*1.5= '+'%.0f'%(x.max()*1.5)
    else:
        xvals=xspace
    ypred = a*(xvals**b)
    ax.plot(xvals,ypred,**kwargs)
    return powfunc

def PowerFit_CI(x,y,xspace='none',ax=plt,**kwargs):
    datadf = pd.DataFrame.from_dict({'x':x,'y':y}).dropna().apply(np.log10) ## put x and y in a dataframe so you can drop ones that don't match up    
    regression = pd.ols(y=datadf['y'],x=datadf['x'])    
    ## Develop power function for x and y
    powfunc = powerfunction(x,y) ## x and y should be Series
    a, b = powfunc['a'].values, powfunc['b'].values
    #print a,b
    if xspace=='none':
        xvals = np.linspace(0,x.max()*1.2)
        #print 'No xspace, calculating xvals: '+str(x.max())+'*1.5= '+str(x.max()*1.5)
    else:
        xvals=xspace
    ypred = a*(xvals**b)
    ax.plot(xvals,ypred,**kwargs)
    ## Confidence interals
    ci=.5
    a_cilo,a_ciup = 10**regression.sm_ols.conf_int(alpha=ci)[1][0],10**regression.sm_ols.conf_int(alpha=ci)[1][1]
    b_cilo,b_ciup = regression.sm_ols.conf_int(alpha=ci)[0][0],regression.sm_ols.conf_int(alpha=ci)[0][1]
    ypred_cilo=a_cilo*(xvals**b_cilo)
    ypred_ciup=a_ciup*(xvals**b_ciup)
    ax.fill_between(xvals,ypred_cilo,ypred_ciup,alpha=0.5,**kwargs)
    return powfunc
## test 
#x= np.linspace(1.0,10.0,10)
#y = 10.0*(x**0.5)
#name = 'x2'
#plt.scatter(x,x2)
#xpowfun = powerfunction(x,x2,name)
#xpower = PowerFit(x,y)
    
def linearfunction(x,y,name='linear rating'):
    datadf = pd.DataFrame.from_dict({'x':x,'y':y}).dropna() ## put x and y in a dataframe so you can drop ones that don't match up    
    datadf = datadf[datadf>=0].dropna() ##verify data is valid (not inf)
    regression = pd.ols(y=datadf['y'],x=datadf['x'],intercept=False)
    pearson = pearson_r(datadf['x'],datadf['y'])[0]
    spearman = spearman_r(datadf['x'],datadf['y'])[0]
    coeffdf = pd.DataFrame({'b':[regression.beta[0]],'r2':[regression.r2],'rmse':[regression.rmse],'pearson':[pearson],'spearman':[spearman]},index=[name])
    return coeffdf
    
def LinearFit(x,y,xspace='none',ax=plt,**kwargs):
    linfunc = linearfunction(x,y)
    a, b = linfunc['a'].values, linfunc['b'].values
    #print a,b
    if xspace=='none':
        xvals = np.linspace(0,x.max()*1.2) ##list of dummy x's as predictors
    else:
        xvals=xspace
    ypred = b*xvals + a ## predicted y from dummy list of x's
    ax.plot(xvals,ypred,**kwargs)
    return linfunc
## test
#x= np.linspace(1.0,10.0,10)
#y = 10*x + 15
#name = 'x2'
#plt.scatter(x,y)
#xlinfun = linearfunction(x,y,name)
#xlinear = LinearFit(x,y)

def nonlinearfunction(x,y,order=2,interceptZero=False):
    datadf = pd.DataFrame.from_dict({'x':x,'y':y}).dropna() ## put x and y in a dataframe so you can drop ones that don't match up    
    datadf = datadf[datadf>=0].dropna() ##verify data is valid (not inf)
    if interceptZero!=True:
        PolyCoeffs = np.polyfit(datadf['x'].values, datadf['y'].values, order) ## calculates polynomial coeffs
        PolyEq = np.poly1d(PolyCoeffs) ## turns the coeffs into an equation
        
    ## Calculate polynomial with a y-intercept of zero    
    if interceptZero==True:
        coeff = np.transpose([datadf['x'].values*datadf['x'].values, datadf['x'].values])
        ((a, b), _, _, _) = np.linalg.lstsq(coeff, datadf['y'].values)
        PolyEq = np.poly1d([a, b, 0])
    return PolyEq

def NonlinearFit(x,y,order=2,interceptZero=False,xspace='none',Ax=plt,**kwargs):
    nonlinfunc = nonlinearfunction(x,y,order,interceptZero)
    #print linfunc
    if xspace=='none':
        xvals = np.linspace(0,x.max()*1.2) ##list of dummy x's as predictors
    else:
        xvals=xspace
    ypred = nonlinfunc(xvals)
    Ax.plot(xvals,ypred,**kwargs)
    return nonlinfunc

## test
#x= np.linspace(1.0,10.0,10)
#y = 10*x + 15
#name = 'x2' 
#plt.scatter(x,x2)
#xnonlinfun = nonlinearfunction(x,x2)
#xnonlinear = NonlinearFit(x,x2)

def plotregressionline(data,ols_object,ax,color):
    slope,intercept = ols_object.beta
    x = np.array([min(data), max(data)])
    y = intercept + slope * x
    ax.plot(x, y,color)
    return
    
def showstormintervals(ax,StormsList,shade_color='grey',show=True):
    ## Storms
    if show==True:
        for storm in StormsList.iterrows(): ## shade over storm intervals
            ax.axvspan(storm[1]['start'],storm[1]['end'],ymin=0,ymax=200,facecolor=shade_color, alpha=0.25)
    return
    
def Sum_Storms(Storm_list,Data,offset=0):
    eventlist = []
    print 'Summing storms...'
    for storm_index,storm in Storm_list.iterrows():
        start = storm['start']-dt.timedelta(minutes=offset) ##if Storms are defined by stream response you have to grab the preceding precip data
        end= storm['end']
        data = True ## Innocent until proven guilty
        try:
            #print str(start) +' '+str(end)
            event = Data.ix[start:end] ### slice list of Data for event
        except KeyError:
            raise
            start = start+dt.timedelta(minutes=15) ## if the start time falls between 2 30minute periods
            print 'change storm start to '+str(start)            
            try:
                event = Data.ix[start:end]
            except KeyError:
                end = end+dt.timedelta(minutes=15)
                print 'change storm end to '+str(end) 
                try:
                    event = Data.ix[start:end]
                except KeyError:
                    print 'no storm data available for storm '+str(start)
                    data = False
                    pass
        if data != False:
            eventcount,eventsum,eventmax=np.nan,np.nan,np.nan
            ## only take events with complete data
            if len(event.dropna()) == len(event):
                eventcount = event.count()
                eventsum = event.sum()
                eventmax = event.max()
                eventlist.append((storm['start'],[storm['start'],storm['end'],eventcount,eventsum,eventmax])) 
            else:
                eventlist.append((storm['start'],[storm['start'],storm['end'],eventcount,eventsum,eventmax])) 
        if data == False:
            eventcount,eventsum,eventmax=np.nan,np.nan,np.nan
            eventlist.append((storm['start'],[storm['start'],storm['end'],eventcount,eventsum,eventmax])) 
    Events=DataFrame.from_items(eventlist,orient='index',columns=['start','end','count','sum','max'])
    return Events
    

def table_to_html_R(dataframe, caption='', table_num='', filename=maindir, save=False, show=True):
    ## convert to R Data Frame
    table_df = com.convert_to_r_dataframe(dataframe)
    ## Send to R
    ro.globalenv['table_df'] = table_df
    ro.globalenv['table_caption'] = 'Table '+str(table_num)+'. '+caption
    ## import htmlTable
    ro.r("library(htmlTable)")
    ro.r("table_out <- htmlTable(table_df, caption=table_caption)")
    htmlcode = com.load_data("table_out")[0]
    if show==True:
        print htmlcode
    if save==True:
        pypandoc.convert(htmlcode, 'html', format='markdown', outputfile= filename)
    return htmlcode



## Year Interval Times
start2012, stop2012 = dt.datetime(2012,1,1,0,0), dt.datetime(2012,12,31,11,59)    
start2013, stop2013 = dt.datetime(2013,1,1,0,0), dt.datetime(2013,12,31,11,59)
start2014, stop2014 = dt.datetime(2014,1,1,0,0), dt.datetime(2015,1,9,11,59)   
#start2015, stop2015 = dt.datetime(2015,1,10,0,0), dt.datetime(2015,4,10,11,59)   
## Field Seasons
fieldstart2012, fieldstop2012 =  dt.datetime(2012,1,5,0,0), dt.datetime(2012,3,29,11,59)    
fieldstart2013, fieldstop2013 =  dt.datetime(2013,2,4,0,0), dt.datetime(2013,7,17,11,59)    
fieldstart2014a, fieldstop2014a =  dt.datetime(2014,1,10,0,0), dt.datetime(2014,3,7,11,59)
fieldstart2014b, fieldstop2014b =  dt.datetime(2014,9,29,0,0), dt.datetime(2015,1,12,11,59)     
## Mitigation
Mitigation = dt.datetime(2014,10,1,0,0)



### LAND COVER
#### Load Land Cover Data
def LandCover_table(browser=True):
    
    ## Read in Excel sheet of data, indexed by subwatershed
    landcover_table = pd.ExcelFile(datadir+'/LandCover/Nuuuli_Watershed_Stats.xlsx').parse('Nuuuli', index_col = 'Watershed')
    landcover_table = landcover_table[['Cumulative Area km2','Cumulative %','Area km2','% of area','% Bare Land','% High Intensity Developed','% Developed Open Space','% Grassland (agriculture)','% Forest','% Scrub/ Shrub','% Disturbed','% Undisturbed']]
    # Format Table data                       
    for column in landcover_table.columns:
        try:
            if column.startswith('%')==True or column == 'Cumulative %':
                #print '1'+ column
                landcover_table.loc[:,column] = landcover_table.loc[:,column]*100.
                landcover_table.loc[:,column] = landcover_table.loc[:,column].round(1)
            else:
                #print '2' + column
                landcover_table.loc[:,column] = landcover_table.loc[:,column].round(2)
        except:
            pass
        
    ## Select the subwatersheds you want
    landcover_table = landcover_table[landcover_table.index.isin(['Nuuuli Upper','Nuuuli Lower','Nuuuli Total'])==True]
    ## Rename the table columns
    landcover_table.columns=['km2 ','% ',' km2',' %','Bare (B)','High Intensity Developed (HI)',
                             'Developed Open Space (DOS)','Grassland (agriculture) (GA)','Forest (F)','Scrub/ Shrub (S)',
                             'Disturbed B+HI+DOS+GA','Undisturbed F+S']
    
    ## convert to R Data Frame
    table_df = com.convert_to_r_dataframe(landcover_table)
    caption="Land use categories in Nu'uuli subwatersheds (NOAA Ocean Service and Coastal Services Center, 2010). Land cover percentages are of the subwatershed."
    table_num=1
    ## Send to R
    ro.globalenv['table_df_vals'] = table_df
    ## format #s
    ro.r("table_df= apply(table_df_vals, 2, function(x) as.character(format(x,digits=2)))")
    ro.r("rownames(table_df)<- rownames(table_df_vals) ")
    #print (ro.r('table_df'))
    ro.globalenv['table_caption'] = 'Table '+str(table_num)+'. '+caption
    ## import htmlTable
    ro.r("library(htmlTable)")
    ## Create table in R
    table_code_str = " \
    table_df, \
    header= c('km<sup>2</sup>','% ','km<sup>2</sup>', '%', '   B   ', '   HI   ', '   DOS   ', '   GA   ', '   F   ', '   S   ', 'Disturbed', 'Undisturbed'), \
    rowlabel='"+landcover_table.index.name+"',\
    align='lccccrrrrrrcc', \
    caption=table_caption, \
    cgroup = c('Cumulative Area','Subwatershed Area','Land cover as % subwatershed area <sup>a</sup>'), \
    n.cgroup = c(2,2,8), \
    tfoot='a. B=Bare, HI=High Intensity Developed, DOS=Developed Open Space, GA=Grassland (agriculture), F=Forest, S=Scrub/Shrub, Disturbed=B+HI+DOS+GA,  Undisturbed=F+S', \
    css.cell = 'padding-left: .5em; padding-right: .2em;'  \
    "
    ## run htmlTable
    ro.r("table_out <- htmlTable("+table_code_str+")")
    ## output to browser
    if browser == True:
        print (ro.r("table_out"))
    ## save to html from R
    ro.r("setwd("+"'"+tabledir+"'"+")")
    ro.r("sink('Table"+str(table_num)+"_Nuuuli_landcover.html')")
    ro.r("print(table_out,type='html',useViewer=FALSE)")
    ro.r("sink()")
    
    ## send back to python
    #htmlcode = com.load_data("table_out")[0] 
    #print htmlcode
    ## output to file through pandoc
    #pypandoc.convert(htmlcode, 'markdown', format='markdown', outputfile= datadir+'landcover.html')
    return landcover_table
#landcover_table = LandCover_table(browser=True)


#### LOAD FIELD DATA
if 'XL' not in locals():
    print 'opening MASTER_DATA excel file...'+dt.datetime.now().strftime('%H:%M:%S')
    XL = pd.ExcelFile(datadir+'MASTER_DATA_NUUULI.xlsx')
if 'XL' in locals():
    print 'MASTER_DATA opened: '+dt.datetime.now().strftime('%H:%M:%S')    


#### Import PRECIP Data
#from precip_data import raingauge#, AddTimu1, AddTimu1Hourly, AddTimu1Daily, AddTimu1Monthly
def raingauge(XL,sheet='',shift=0):
    print 'loading precip: '+sheet+'...'
    #my_parser= lambda x: dt.datetime.strptime(x,"%m/%d/%Y %H:%M")
    gauge = XL.parse(sheet,header=1,index_col=0,parse_cols='B,C',parse_dates=True)#,date_parser=my_parser)
    gauge= gauge.shift(shift)
    gauge = gauge*0.254 ##hundredths to mm
    gauge.columns=['Events']
    return gauge    


if 'Precip' not in locals():
    ## Timu-Fagaalu 1 (by the Quarry)
    ## 2012-2013 Data
    Precip = raingauge(XL,'Timu-Nuuuli1',0) ## (path,sheet,shift) no header needed    
    #Precip= Precip.reindex(pd.date_range(dt.datetime(2012,1,6,17,51),dt.datetime(2013,12,31,23,59),freq='1Min'))
    Precip.columns=['Timu1']
    Precip['Timu1-15']=Precip['Timu1'].resample('15Min',how='sum')
    Precip['Timu1-30']=Precip['Timu1'].resample('30Min',how='sum')
    # Hourly
    Precip['Timu1hourly']= Precip['Timu1'].resample('H',how='sum')
    Precip['Timu1hourly'].dropna().to_csv(datadir+'OUTPUT/Timu1hourly.csv',header=['Timu1hourly'])
    # Daily
    Precip['Timu1daily'] = Precip['Timu1'].resample('D',how='sum')
    Precip['Timu1daily'].dropna().to_csv(datadir+'OUTPUT/Timu1daily.csv',header=['Timu1daily'])
    # Monthly
    Precip['Timu1monthly'] = Precip['Timu1'].resample('MS',how='sum') ## Monthly Precip
    Precip['Timu1monthly'].dropna().to_csv(datadir+'OUTPUT/Timu1monthly.csv',header=['Timu1monthly'])
    
    ## Timu-Fagaalu2 
    Precip['Timu2']=raingauge(XL,'Timu-Nuuuli2',0)
    Precip['Timu2-15']=Precip['Timu2'].resample('15Min',how='sum')
    Precip['Timu2-30']=Precip['Timu2'].resample('30Min',how='sum')
    # Hourly
    Precip['Timu2hourly']= Precip['Timu2-30'].resample('H',how='sum')
    Precip['Timu2hourly'].dropna().to_csv(datadir+'OUTPUT/Timu2hourly.csv',header=['Timu2hourly'])
    # Daily
    Precip['Timu2daily'] = Precip['Timu2-30'].resample('D',how='sum')
    Precip['Timu2daily'].dropna().to_csv(datadir+'OUTPUT/Timu2daily.csv',header=['Timu2daily'])
    # Monthly
    Precip['Timu2monthly'] = Precip['Timu2-30'].resample('MS',how='sum') ## Monthly Precip
    Precip['Timu2monthly'].dropna().to_csv(datadir+'OUTPUT/Timu2monthly.csv',header=['Timu2monthly'])


#### Import BAROMETRIC Data: Fagaalu Weather Station

#from load_from_MASTER_XL import WeatherStation
def WeatherStation(XL,sheet=''):
    print 'loading Wx: '+sheet+'...'
    ## Fagaalu Weather Station
    #my_parser= lambda x,y: dt.datetime.strptime(x+y,"%m/%d/%Y%I:%M %p")
    Wx= XL.parse(sheet,skiprows=1,header=0,parse_cols='A:AD',parse_dates=[['Date','Time']],index_col=['Date_Time'],na_values=['---'])
    Wx.columns=['TempOut', 'HiTemp', 'LowTemp', 'OutHum', 'DewPt', 'WindSpeed', 'WindDir', 'WindRun', 'HiSpeed', 'HiDir', 'WindChill', 'HeatIndex', 'THWIndex', 'Bar', 'Rain', 'RainRate', 'HeatD-D', 'CoolD-D', 'InTemp', 'InHum', 'InDew', 'InHeat', 'InEMC', 'InAirDensity', 'WindSamp', 'WindTx', 'ISSRecept', 'Arc.Int.']
    return Wx
    
if 'FP' not in locals():
    print 'Opening Weather Station data...'
    FP = WeatherStation(XL,'FP-15min')
    



#### Import BAROMETRIC Data: NDBC
##load data from NDBC NSTP6 station at DMWR, Pago Harbor
## To get more NSTP6 data either go to their website and copy and paste the historical data
## or use wundergrabber_NSTP6-REALTIME.py and copy and paste frome the .csv
def ndbc(datafile = datadir+'BARO/NSTP6/NSTP6-2012-14.xlsx'):
    print 'Loading NDBC NSTP6 barometric data...'
    try:
        ndbc_data = pd.DataFrame.from_csv(datadir+'BARO/NSTP6/NDBC_Baro.csv')
    except:
        ndbcXL = pd.ExcelFile(datafile)
        ndbc_parse = lambda yr,mo,dy,hr,mn: dt.datetime(yr,mo,dy,hr,mn)
        ndbc_data = ndbcXL.parse('NSTP6-2012-14',header=0,skiprows=1,parse_dates=[['#yr','mo','dy','hr','mn']],index_col=0,date_parser=ndbc_parse,
                             na_values=['9999','999','99','99.0'])
        ndbc_data.to_csv(datadir+'BARO/NSTP6/NDBC_Baro.csv')
    #local = pytz.timezone('US/Samoa')
    #ndbc_data.index = ndbc_data.index.tz_localize(pytz.utc).tz_convert(local)
    print 'NDBC loaded'
    return ndbc_data

NDBCbaro = ndbc(datafile = datadir+'BARO/NSTP6/NSTP6-2012-14.xlsx')
NDBCbaro = NDBCbaro['hPa'].resample('15Min')
NDBCbaro = NDBCbaro.interpolate(method='linear',limit=4)
NDBCbaro.columns=['NDBCbaro']
NDBCbaro=NDBCbaro.shift(-44) ## UTC to Samoa local  =11 hours =44x15min
NDBCbaro = NDBCbaro-.022

## Barologger 
def barologger(XL,sheet=''):
    print 'loading Wx: '+sheet+'...'
    ## Fagaalu Weather Station
    #my_parser= lambda x,y: dt.datetime.strptime(x+y,"%m/%d/%Y%I:%M %p")
    Baro=  XL.parse(sheet,header=11,parse_cols='A,B,D',parse_dates=[['Date','Time']],index_col=['Date_Time'])
    Baro.columns= ['LEVEL']
    Baro=Baro.resample('15Min',how='mean')
    return Baro
    

 
## Build data frame of barometric data: Make column 'baropress' with best available data
 ## Fill priority = FP,NDBC,TAFUNA,TULA (TAFUNA and TULA have been deprecated, reside in other scripts)
allbaro = pd.DataFrame(NDBCbaro/10).reindex(pd.date_range(start2012,stop2014,freq='15Min'))
allbaro['FPbaro']=FP['Bar']/10
allbaro['NDBCbaro']=NDBCbaro/10
#allbaro['BaroLogger']=BaroLogger

## create a new column and fill with FP or NDBC
allbaro['Baropress']=allbaro['FPbaro'].where(allbaro['FPbaro']>0,allbaro['NDBCbaro']) 



#### Import PT Data
# ex. PT_Levelogger(allbaro,PTname,datapath,tshift=0,zshift=0): 
#from load_from_MASTER_XL import PT_Hobo,PT_Levelogger
def PT_Hobo(allbaro,PTname,XL,sheet='',tshift=0,zshift=0): # tshift in 15Min(or whatever the timestep is), zshift in cm
    print 'loading HOBO PT: '+sheet+'...'
    PT = XL.parse(sheet,header=1,index_col=0,parse_cols='B,C',parse_dates=True)
    PT.columns=['Pressure']
    PT=PT.resample('15Min',how='mean')
    PT=PT.shift(tshift) ## shift by 3 hours (12 x 15minutes)
    PT['barodata']=allbaro['Baropress']
    PT['stage(cm)']=(PT['Pressure']-PT['barodata'])*.102*100.0 ## hPa  to cm
    #PT['stage']=PT['stage'].where(PT['stage']>0,PT['barodata']) ## filter negative values
    PT['stage(cm)']=PT['stage(cm)'].round(1)  
    PT['stage(cm)']=PT['stage(cm)']+zshift
    PT['Uncorrected_stage']=PT['stage(cm)'].round(0)
    return PT

def PT_Levelogger(allbaro,PTname,XL,sheet,tshift=0,zshift=0): # tshift in hours, zshift in cm
    print 'loading Levelogger PT: '+sheet+'...'
    PT = XL.parse(sheet,header=11,parse_cols='A,B,D',parse_dates=[['Date','Time']],index_col=['Date_Time'])
    PT.columns= ['LEVEL']
    PT=PT.resample('15Min',how='mean')
    PT['barodata']=allbaro['Baropress']
    PT=PT.shift(tshift) ## shift by 3 hours (12 x 15minutes)
    PT['stage(cm)']=(PT['LEVEL']-PT['barodata'])*.102*100.0
    #PT['stage']=PT['stage'].where(PT['stage']>0,0) ## filter negative values
    PT['stage(cm)']=PT['stage(cm)'].round(1)  
    PT['stage(cm)']=PT['stage(cm)']+zshift
    PT['Uncorrected_stage']=PT['stage(cm)'].round(0)
    return PT


## LOAD PT data
# PT1 = Downstream
PT1 = PT_Levelogger(allbaro,'PT-Nuuuli1',XL,'PT-Nuuuli1',0)
PT1['stage(cm)'][PT1['stage(cm)']<=0] = np.nan

# PT2 = Downstream
PT2 = PT_Levelogger(allbaro,'PT-Nuuuli2',XL,'PT-Nuuuli2',0)
PT2['stage(cm)'][PT2['stage(cm)']<=0] = np.nan


## Reindex
PT1 = PT1.reindex(pd.date_range(start2013,stop2014,freq='15Min'))
PT2 = PT2.reindex(pd.date_range(start2013,stop2014,freq='15Min'))

def plot_uncorrected_stage_data(show=False):
    fig, (baro,pt1,pt2,t) = plt.subplots(4,sharex=True,sharey=False,figsize=(12,6))
    for ax in [baro,pt1,pt2]:
        ax.xaxis.set_visible(False)
    allbaro['Baropress'].plot(ax=baro,c='k',label='Barometric Pressure (kPa)')
    #allbaro['Baropress'].plot(ax=t,c='k',label='Barometric Pressure (kPa)')
    baro.set_ylabel('kPa')
    baro.legend()
    ## PT1 Downstream
    PT1list = [PT1]
    count = 0
    count_dict = {1:'aa',2:'ab',3:'ba',4:'bb',5:'bc',6:'c'}
    count_dict = {1:' '}
    for PT in PT1list:
        count+=1
        try:
            PT['Pressure'].plot(ax=pt1,label='PT1'+count_dict[count],c='r')
            PT['stage(cm)'].plot(ax=t,label='PT1'+count_dict[count],c='r')
        except KeyError:
            PT['LEVEL'].plot(ax=pt1,label='PT1'+count_dict[count],c='r')
            PT['stage(cm)'].plot(ax=t,label='PT1'+count_dict[count],c='r')
    pt1.legend()
    pt1.set_ylabel('kPa')
    ## PT2 Upstream
    PT2list = [PT2]
    count = 0
    count_dict = {1:'aa',2:'ab',3:'b',4:'c',5:'d',6:'e',7:'f',8:'g'}
    count_dict = {1:' '}
    for PT in PT2list:
        count+=1
        try:
            PT['Pressure'].plot(ax=pt2,label='PT2'+count_dict[count],c='g')
            PT['stage(cm)'].plot(ax=t,label='PT2'+count_dict[count],c='g')
        except KeyError:
            PT['LEVEL'].plot(ax=pt2,label='PT2'+count_dict[count],c='g')
            PT['stage(cm)'].plot(ax=t,label='PT2'+count_dict[count],c='g')
    pt2.legend()
    pt2.set_ylabel('kPa')
    t.legend(ncol=3)
    t.set_ylabel('cm')
    
    if show==True:
        plt.tight_layout(pad=0.01)
        plt.show()
    return
#plot_uncorrected_stage_data(show=True)

## STAGE DATA FOR PT's
#### FINAL STAGE DATA with CORRECTIONS
Nuuuli_stage_data = pd.DataFrame({'N1':PT1['stage(cm)'],'N2':PT2['stage(cm)']})
Nuuuli_stage_data = Nuuuli_stage_data.reindex(pd.date_range(start2013,stop2014,freq='15Min'))


#### STAGE TO DISCHARGE ####
#from stage2discharge_ratingcurve import AV_RatingCurve#, calcQ, Mannings_rect, Weir_rect, Weir_vnotch, Flume


### Calculate Q from a single AV measurement
#fileQ = calcQ(datadir+'Q/LBJ_4-18-13.txt','LBJ',Fagaalu_stage_data,slope=Slope,Mannings_n=n,trapezoid=True)
## and save to CSV
#pd.concat(fileQ).to_csv(datadir+'Q/LBJ_4-18-13.csv')

### Area Velocity and Mannings from in situ measurments
## Returns DataFrame of Stage (cm) and Discharge (L/sec) calc. from AV measurements with time index

def Stage_Q_AV_RatingCurve(path,location,stage_data,slope=.01,Mannings_n=.033,trapezoid=True,printResults=False):
    Filelist = os.listdir(path)
    ## iterate over files in directory to get Flow.txt file
    for f in Filelist:
        ## Select Flow.txt file
        if f.endswith('Flow.txt')==True and f.startswith(location)==True:
            print 'AV measurements file selected for analysis: '+f
            ## Open File, create blank parameters
            Flowfile = open(path+f)
            Qdf = pd.DataFrame() ## empty dataframe to append calculated Q
            for line in Flowfile:
                split = line.strip('\n').split('\t')
                #print split
                # Test if data is number
                try:
                    a= float(split[0]) ## dummy test
                    isfloat=True
                except ValueError:
                    isfloat=False            
                ## Determine DateTime of AV measurment
                if split[0]==location:
                    ## Create empty dataframe
                    df = pd.DataFrame(columns=['dist','depth','flow']) ## empty dataframe for Flowmeter data
                    date, time = split[1].split('/'),split[2]
                    if len(time)==3:
                        time = '0'+time
                    DateTime = dt.datetime(int(date[2]),int(date[0]),int(date[1]),int(time[0:2]),int(time[2:]))
                    DateTime = RoundTo15(DateTime)
                    #print DateTime
                ## Append data
                elif isfloat==True:
                    df=df.append(pd.DataFrame({'dist':split[0],'depth':split[1],'flow':split[2]},index=[DateTime]))
                elif split[0]=='Location' or split[0]=='Field Measurements' or split[0]=='Dist(S to N)(ft)' or split[0]=='Dist(ft)':
                    pass
                
                ## At the end of that AV measurment, calculate Q
                elif split[0]=='-':
                    #print 'calculating Q for '+str(DateTime)
                    df = df.astype('float')
                    if trapezoid==True:
                        ## Depth/flow measurements are made at the midpoint of the trapezoid, dimensions of the trapezoid have to be determined
                        df['right']= df['dist'].shift(-1).sub(df['dist'],fill_value=0) ##Distance next - Distance = the width to the right
                        df['left'] = df['dist'].sub(df['dist'].shift(1),fill_value=0) ## Distance previous - Distance = the width to the left
                        df['right'] = df['right'].where(df['right']>0,0)
                        df['width']=((df['right'] + df['left'])/2)*12*2.54 ## 2nd mark - first; then convert to cm
                        df['b1']=(df['depth'].add(df['depth'].shift(1),fill_value=0))/2 ## gives average of Depth Above and depth
                        df['b2']=(df['depth'].add(df['depth'].shift(-1),fill_value=0))/2 ## gives average of Depth Below and depth
                        ## Formula for area of a trapezoid = 1/2 * (B1+B2) * h; h is width of the interval and B1 and B2 are the depths at the midpoints between depth/flow measurements
                        df['trapezoidal-area-cm2']=.5*(df['b1']+df['b2'])*df['width'] 
                        df['trapezoidal-area-m2']=df['trapezoidal-area-cm2']/10000 ##cm2 to m2
                        df['AV']=df['trapezoidal-area-m2']*df['flow'] *1000 ## m2 x m/sec x 1000 = L/sec
                        AV_Q = df['AV'].sum()
                        Area = df['trapezoidal-area-m2'].sum()
                        V = df['flow'].mean()
                        ## Wetted perimeter doesn't use midpoints between depth/flow measurments
                        ## WP = SQRT((Dnext-D)^2 + Width^2)
                        df['WP']=  ((df['depth'].sub(df['depth'].shift(1),fill_value=0))**2 + (df['dist'].sub(df['dist'].shift(1),fill_value=0)*12*2.54)**2)**0.5 
                        ## cm to m; and only take WP values where the depth to the left is not zero
                        df['WP']=(df['WP']*(df['b1']>0))/100 
                        
                        WP = df['WP'].sum()
                        R = Area / WP ## m2/m = m
                        ## Mannings = 1/n * R^2/3 * S^1/2
                        S = slope
                        ## Jarrett (1990) equation for n
                        ## n = 0.32*(S**0.30)*(R**-0.16)
                        if Mannings_n == 'Jarrett':
                            n = 0.32*(S**0.30)*(R**-0.16)
                        else:
                            n = Mannings_n
                        ManningV = (1.0/n) * (R**(2.0/3.0)) * (S**0.5)
                        ManningQ = ManningV * Area * 1000 ## L/Sec
                    elif trapezoid==False:
                        df = df.set_value(len(df),'dist',df['dist'][-1]) ## add a dummy distance value
                        valbelow = df['dist'].shift(-1).sub(df['dist'],fill_value=0) ## Width is value below - dist value
                        valabove = df['dist'].sub(df['dist'].shift(1),fill_value=0)
                        df['width']=(valbelow.add(valabove)/2)*12*2.54 ## 2nd mark - first
                        df['rectangular-area']=df['depth']*(df['width'])/10000 ##cm2 to m2
                        df['AV']=df['rectangular-area']*df['flow']
                        
                        AV_Q = df['AV'].sum().round(0)
                        Area = df['rectangular-area'].sum()
                        V = df['flow'].mean()
                        ManningV = np.nan
                    try:
                        stage = stage_data[location].ix[DateTime] ## Get Stage data
                        print location+' '+str(DateTime)+' '+'Stage= '+str(stage)+' Q= '+str(AV_Q)
                    except:
                        stage =np.nan
                        print location+' '+str(DateTime)+' '+'Stage= '+str(stage)+' Q= '+str(AV_Q)
                    Qdf = Qdf.append(pd.DataFrame({'stage(cm)':stage,'Q-AV(L/sec)':round(AV_Q,0),'Q-AManningV(L/sec)':round(ManningQ,0),
                    'Area(m2)':Area,'V(m/s)':V,'ManningV(m/s)':ManningV,'WP':WP,'R':R,'n':Mannings_n},index=[DateTime]))
                    print str(DateTime)+': stage='+'%.2f'%stage+' Q= '+'%.0f'%AV_Q+' ManningQ(L/sc)= '+'%.2f'%ManningQ+' Area(m2)='+'%.2f'%Area+' WP='+'%.2f'%WP+' R='+'%.2f'%R+' n='+'%.3f'%Mannings_n
                    if printResults==True:                    
                        print df              
    return Qdf[['stage(cm)','Q-AV(L/sec)','Q-AManningV(L/sec)','V(m/s)','ManningV(m/s)','Area(m2)','R','WP','n']]


### Discharge using Mannings and Surveyed Cros-section
#from ManningsRatingCurve import Mannings, Mannings_Series
def Mannings_Q_from_CrossSection(Cross_section_file,sheetname,Slope,Manning_n,k=1,stage_start=.01,stage_end=None,show=False,save=False,filename=''):    
    ## Open and parse file; drop NA  
    print Cross_section_file+' '+sheetname
    print 'Slope: '+str(Slope)+' Mannings n: '+str(Manning_n)
    XL = pd.ExcelFile(Cross_section_file) 
    df = XL.parse(sheetname,header=4,parse_cols='F:H')
    df = df.dropna()
    ## Mannings Parameters S:slope, n:Mannings n
    S = Slope # m/m
    n= Manning_n
    ## empty lists
    areas, wp, r, Man_n, v, q, = [],[],[],[],[],[]
    ## Stage data
    ## one stage measurement
    if stage_end == None:
        print 'Stage: '+str(stage_start)
        stages = np.array([stage_start])
    ## start and end stage
    elif stage_start != stage_end:
        print 'Stage_start: '+str(stage_start)+' Stage_end: '+str(stage_end)
        stages = np.arange(stage_start,stage_end,.1) #m
    ## stage Series         
    elif type(stage_start)==pd.Series:
        print 'Stage Series...'
        stages = stage_start.to_list()
        
    for stage in stages:
        print 'stage: '+str(stage)
        df['y1'] = df['depth']+df['Rod Reading'].max()
        df['y2'] = stage
        df['z'] = df['y2']-df['y1']
        df['z'] = df['z'][df['z']>=0]
        
        x = df['Dist'].values
        y1 = df['y1'].values
        y2 = df['y2'].values
        
        z = y2-y1
        z= np.where(z>=0,z,0)
        Area = np.trapz(z,x)
        
        ## Wetted Perimeter
        df['dx'] = df['Dist'].sub(df['Dist'].shift(1),fill_value=0)
        df['dy'] = df['z'].sub(df['z'].shift(1),fill_value=0)
        df['wp'] = (df['dx']**2 + df['dy']**2)**0.5
        print df        
        
        WP = df['wp'].sum()
        R = (Area/WP) ## m2/m = m
        ## Jarrett (1990) equation for n
        ## n = 0.32*(S**0.30)*(R**-0.16)
        if Manning_n == 'Jarrett':
            n = 0.32*(S**0.30)*(R**-0.16)
            n= n *k
        ## Mannings = (1/n * R^2/3 * S^1/2)
        ManningV = 1/n * (R**(2.0/3.0)) * S**0.5
        ManningQ = ManningV * Area ## M3/s
        
        plt.ioff()          
        fig, ax1 = plt.subplots(1)
        ax1.plot(df['Dist'],df['y1'],'-o',c='k')
        ax1.fill_between(df['Dist'], df['y1'], stage,where = df['y1']<=stage,alpha=.5, interpolate=True)
        
        ax1.annotate('stage: '+'%.2f'%stage+'m',xy=(0,1.5+.45))
        ax1.annotate('Mannings n: '+'%.3f'%n,xy=(0,1.5+.03))
        ax1.annotate('Area: '+'%.3f'%Area+'m2',xy=(0,1.5+.25))
        ax1.annotate('WP: '+'%.2f'%WP+'m',xy=(df['Dist'].mean(),1.5+.03))
        ax1.annotate('Manning V: '+'%.2f'%ManningV+'m/s ',xy=(df['Dist'].mean(),1.5+.25))
        ax1.annotate('Manning Q: '+'%.3f'%ManningQ+'m3/s',xy=(df['Dist'].mean(),1.5+.45))
        plt.axes().set_aspect('equal')
        plt.xlim(-1,df['Dist'].max()+1),plt.ylim(-1,2 + 1.)
    
        areas.append(Area)
        wp.append(WP)
        r.append(R)
        Man_n.append(n)
        v.append(ManningV)
        q.append(ManningQ)
        show_plot(show,fig)
        savefig(save,filename) 
        plt.close('all')
        plt.ion()
    
    DF = pd.DataFrame({'stage(m)':stages,'area(m2)':areas,'wp(m)':wp,'r':r,'Man_n':Man_n,'vel(m/s)':v,'Q(m3/s)':q}) 
    return DF,df

def Mannings_Q_from_stage_data(Cross_section_file,sheetname,stage_data,Slope,Manning_n,k=1):    
    ## Open and parse file; drop NA  
    print Cross_section_file+' '+sheetname
    print 'Slope: '+str(Slope)+' Mannings n: '+str(Manning_n)
    XL = pd.ExcelFile(Cross_section_file) 
    df = XL.parse(sheetname,header=4,parse_cols='F:H')
    df = df.dropna()
    ## Mannings Parameters S:slope, n:Mannings n
    S = Slope # m/m
    n= Manning_n
    ## empty lists
    areas, wp, r, Man_n, v, q, = [],[],[],[],[],[]
    ## Stage data
    stage_data = stage_data/100 ## cm to m
    for stage in stage_data.values:
        print 'stage: '+str(stage)
        df['y1'] = df['depth']+df['Rod Reading'].max()
        df['y2'] = stage
        df['z'] = df['y2']-df['y1']
        df['z'] = df['z'][df['z']>=0]
        x = df['Dist'].values
        y1 = df['y1'].values
        y2 = df['y2'].values
        z = y2-y1
        z= np.where(z>=0,z,0)
        Area = np.trapz(z,x)
        ## Wetted Perimeter0.01
        df['dx'] =df['Dist'].sub(df['Dist'].shift(1),fill_value=0)
        df['dy'] = df['z'].sub(df['z'].shift(1),fill_value=0)
        df['wp'] = (df['dx']**2 + df['dy']**2)**0.5
        WP = df['wp'].sum()
        R = (Area/WP) ## m2/m = m
        ## Jarrett (1990) equation for n
        ## n = 0.32*(S**0.30)*(R**-0.16)
        if Manning_n == 'Jarrett':
            n = 0.32*(S**0.30)*(R**-0.16) 
            n = n * k
        ## Mannings = (1R^2/3 * S^1/2)/n
        ManningV = 1/n * (R**(2.0/3.0)) * (S**0.5) ## m/s
        ManningQ = ManningV * Area ## M3/s
        ManningQ= round(ManningQ,3)
        areas.append(Area)
        wp.append(WP)
        r.append(R)
        Man_n.append(n)
        v.append(ManningV)
        q.append(ManningQ)        
    DF = pd.DataFrame({'stage(m)':stage_data.values,'area(m2)':areas,'wp(m)':wp,'r':r,'Man_n':Man_n,'vel(m/s)':v,'Q(m3/s)':q},index=stage_data.index)
    return DF

## N1 AV measurements
## Mannings parameters for A-ManningV
N1_Slope = 0.0013 # m/m
N1_n= .067 # Mountain stream rocky bed and rivers with variable sections and veg along banks (Dunne 1978)
N1_k = 1
#DataFrame with Q from AV measurements, Q from measured A with Manning-predicted V, stage, and Q from Manning's and assumed rectangular channel
N1stageDischarge = Stage_Q_AV_RatingCurve(datadir+'Q/Flow_Files/','N1',Nuuuli_stage_data,slope=N1_Slope,Mannings_n=N1_n,trapezoid=True)#.dropna() 
N1stageDischargeLog = N1stageDischarge.apply(np.log10) #log-transformed version

## N1: Discharge Ratings
## Linear
N1_AV= pd.ols(y=N1stageDischarge['Q-AV(L/sec)'],x=N1stageDischarge['stage(cm)'],intercept=True) 
## Power
N1_AVLog= pd.ols(y=N1stageDischargeLog['Q-AV(L/sec)'],x=N1stageDischargeLog['stage(cm)'],intercept=True) #linear fit to log-transformed stage and Q
## Linear with Mannings and measured Area
N1_AManningV = pd.ols(y=N1stageDischarge['Q-AManningV(L/sec)'],x=N1stageDischarge['stage(cm)'],intercept=True)
## Power with Mannings and measured Area
N1_AManningVLog = pd.ols(y=N1stageDischargeLog['Q-AManningV(L/sec)'],x=N1stageDischargeLog['stage(cm)'],intercept=True)

## Mannings from Cross-Section
## Max/Min stage
N1_max_stage, N1_min_stage = PT1['stage(cm)'].max(), PT1['stage(cm)'].min()
## Max
#Mannings_Q_from_CrossSection(datadir+'Q/Cross_Section_Surveys/N1_cross_section.xlsx','N1-2_m',Slope=N1_Slope,Manning_n=N1_n,k=N1_k,stage_start=N1_max_stage/100,show=True,save=True,filename=datadir+'Q/Cross_Section_Surveys/N1_Mannings_Cross_Section_max_stage')
## Min
#Mannings_Q_from_CrossSection(datadir+'Q/Cross_Section_Surveys/N1_cross_section.xlsx','N1-2_m',Slope=N1_Slope,Manning_n=N1_n,k=N1_k,stage_start=N1_min_stage/100,show=True,save=True,filename=datadir+'Q/Cross_Section_Surveys/N1_Mannings_Cross_Section_min_stage')

## Read N1_Man Discharge from .csv, or calculate new if needed
if 'N1_Man' not in locals():
    try:
        print 'Loading Mannings Q for DAM from CSV'
        N1_Man_reduced = pd.DataFrame.from_csv(datadir+'Q/Manning_Q_files/N1_Man_reduced.csv')
        N1_Man = pd.DataFrame.from_csv(datadir+'Q/Manning_Q_files/N1_Man.csv')
    except:
        print 'Calculate Mannings Q for N1 and saving to CSV'
        N1_stage_reduced = Nuuuli_stage_data['N1'].dropna().round(0).drop_duplicates().order()
        N1_stage_reduced = N1_stage_reduced[N1_stage_reduced>0]
        N1_Man_reduced = Mannings_Q_from_stage_data(datadir+'Q/Cross_Section_Surveys/N1_cross_section.xlsx','N1-2_m',Slope=N1_Slope,Manning_n=N1_n,k=N1_k,stage_data=N1_stage_reduced)
        N1_Man_reduced.to_csv(datadir+'Q/Manning_Q_files/N1_Man_reduced.csv')
        N1_stage= Nuuuli_stage_data['N1']
        N1_Man= Mannings_Q_from_stage_data(datadir+'Q/Cross_Section_Surveys/N1_cross_section.xlsx','N1-2_m',Slope=N1_Slope,Manning_n=N1_n,k=N1_k,stage_data=N1_stage)
        N1_Man.to_csv(datadir+'Q/Manning_Q_files/N1_Man.csv')
        pass

def Manning_AV_r2(ManningsQ_Series,AV_Series):
    # LBJ Mannings = y predicted
    ManQ, Manstage = ManningsQ_Series['Q(m3/s)']*1000., ManningsQ_Series['stage(m)']*100.
    y_predicted = pd.DataFrame({'Q_Man':ManQ.values},index=Manstage).sort()
    y_predicted.index.name = 'stage(cm)'
    ## LBJ AV  = y
    AV_Q, AVstage = AV_Series['Q-AV(L/sec)'], AV_Series['stage(cm)'].apply(np.int)
    y_ = pd.DataFrame({'Q_AV':AV_Q.values}, index=AVstage.values).sort()
    y_['idx'] = y_.index
    y_ = y_[y_.duplicated(cols='idx') == False]
    y_['Q_Man'] = y_predicted['Q_Man']
    y_=  y_.dropna() # keep it clean
    
    #y_.loc[28,'Q_Man'] = 172
    y_bar = y_['Q_AV'].mean()
    y_var = (y_['Q_AV'] - y_bar)**2.
    ss_tot = y_var.sum()
    y_res = (y_['Q_AV']-y_['Q_Man'])**2.
    ss_res = y_res.sum()
    r2 = 1-(ss_res/ss_tot)
    return  r2
#N1_Man_r2 = Manning_AV_r2(N1_Man_reduced,N1stageDischarge)
N1_Man_r2 = 'N/A'

    

### Compare Discharge Ratings from different methods
def plotQratingN1(ms=6,show=False,log=False,save=False,filename=figdir+''): ## Rating Curves
    mpl.rc('lines',markersize=ms)
    title="Water Discharge Ratings for N1"
    fig, (site_N1, site_N1_zoom)= plt.subplots(1,2,figsize=(8,4))
    xy = np.linspace(0,8000,8000)
    
    site_N1.text(0.1,0.95,'(a)',verticalalignment='top', horizontalalignment='right',transform=site_N1.transAxes,color='k',fontsize=10,fontweight='bold')
    site_N1_zoom.text(0.1,0.95,'(b)',verticalalignment='top', horizontalalignment='right',transform=site_N1_zoom.transAxes,color='k',fontsize=10,fontweight='bold')
    
    #N1 AV Measurements and Rating Curve
    site_N1.plot(N1stageDischarge['Q-AV(L/sec)'][start2013:stop2013],N1stageDischarge['stage(cm)'][start2013:stop2013],'^',color='k',fillstyle='none',label='AV 2013') 
    site_N1.plot(N1stageDischarge['Q-AV(L/sec)'][start2014:stop2014],N1stageDischarge['stage(cm)'][start2014:stop2014],'s',color='k',fillstyle='none',label='AV 2014') 
    site_N1_zoom.plot(N1stageDischarge['Q-AV(L/sec)'][start2013:stop2013],N1stageDischarge['stage(cm)'][start2013:stop2013],'^',color='k',fillstyle='none',label='AV 2013') 
    site_N1_zoom.plot(N1stageDischarge['Q-AV(L/sec)'][start2014:stop2014],N1stageDischarge['stage(cm)'][start2014:stop2014],'s',color='k',fillstyle='none',label='AV 2014') 

    ## N1 MODELS
    ## Mannings for AV measurements
    site_N1.plot(N1stageDischarge['Q-AManningV(L/sec)'],N1stageDischarge['stage(cm)'],'.',ls='None',c='grey',label='A-ManningsV')
    site_N1_zoom.plot(N1stageDischarge['Q-AManningV(L/sec)'],N1stageDischarge['stage(cm)'],'.',ls='None',c='grey',label='A-ManningsV')
    ## N1 Power
    N1_AVpower = powerfunction(N1stageDischarge['Q-AV(L/sec)'],N1stageDischarge['stage(cm)'])    
    PowerFit(N1stageDischarge['Q-AV(L/sec)'],N1stageDischarge['stage(cm)'],xy,site_N1,c='grey',ls='-',label='AV power law '+r'$r^2$'+"%.2f"%N1_AVpower['r2'])    
    PowerFit(N1stageDischarge['Q-AV(L/sec)'],N1stageDischarge['stage(cm)'],xy,site_N1_zoom,c='grey',ls='-',label='AV power law '+r'$r^2$'+"%.2f"%N1_AVpower['r2'])        
    ## N1 Mannings from stream survey
    N1_ManQ, N1_Manstage = N1_Man_reduced['Q(m3/s)']*1000, N1_Man_reduced['stage(m)']*100
    site_N1.plot(N1_ManQ, N1_Manstage,'-',markersize=2,c='k',label='Mannings: n='+str(N1_n)+r'$ r^2$'+N1_Man_r2)
    site_N1_zoom.plot(N1_ManQ,N1_Manstage,'-',markersize=2,c='k',label='Mannings')
    ## Label point -click
    labelindex_subplot(site_N1, N1stageDischarge.index,N1stageDischarge['Q-AV(L/sec)'],N1stageDischarge['stage(cm)'])
    labelindex_subplot(site_N1_zoom, N1stageDischarge.index,N1stageDischarge['Q-AV(L/sec)'],N1stageDischarge['stage(cm)'])
    ## Label subplots    
    site_N1.set_ylabel('Stage(cm)'),site_N1.set_xlabel('Q(L/sec)'),site_N1_zoom.set_xlabel('Q(L/sec)')
    ## Format subplots
    site_N1.set_ylim(0,PT1['stage(cm)'].max()+10)#,site_N1.set_xlim(0,N1_AVnonLinear(PT1['stage'].max()+10))
    site_N1_zoom.set_ylim(0,45), site_N1_zoom.set_xlim(0,1600)
    ## Legends
    site_N1.legend(loc='lower right',fancybox=True)  
    ## Figure title
    #plt.suptitle(title,fontsize=16)
    fig.canvas.manager.set_window_title('Figure : '+title) 
    logaxes(log,fig)
    for ax in fig.axes:
        ax.autoscale_view(True,True,True)
    plt.tight_layout(pad=0.1)
    show_plot(show,fig)
    savefig(save,filename)
    return
plotQratingN1(ms=6,show=True,log=False,save=True,filename=figdir+'Stage-Q rating N1')
#plotQratingN1(ms=6,show=True,log=True,save=False,filename=figdir+'')

## N2 AV measurements
## Mannings parameters for A-ManningV
N2_Slope = 0.0023 # m/m
N2_n = .067 # Mountain stream rocky bed and rivers with variable sections and veg along banks (Dunne 1978)
N2_k = 1
#DataFrame with Q from AV measurements, Q from measured A with Manning-predicted V, stage, and Q from Manning's and assumed rectangular channel
N2stageDischarge = Stage_Q_AV_RatingCurve(datadir+'Q/Flow_Files/','N2',Nuuuli_stage_data,slope=N2_Slope,Mannings_n=N2_n,trapezoid=True).dropna() 
N2stageDischargeLog = N2stageDischarge.apply(np.log10) #log-transformed version


## N2: Discharge Ratings
## Linear
N2_AV= pd.ols(y=N2stageDischarge['Q-AV(L/sec)'],x=N2stageDischarge['stage(cm)'],intercept=False) 
## Power
N2_AVLog= pd.ols(y=N2stageDischargeLog['Q-AV(L/sec)'],x=N2stageDischargeLog['stage(cm)'],intercept=False) #linear fit to log-transformed stage and Q
## Linear with Mannings and measured Area
N2_AManningV = pd.ols(y=N2stageDischarge['Q-AManningV(L/sec)'],x=N2stageDischarge['stage(cm)'],intercept=False)
## Power with Mannings and measured Area
N2_AManningVLog = pd.ols(y=N2stageDischargeLog['Q-AManningV(L/sec)'],x=N2stageDischargeLog['stage(cm)'],intercept=False)

## Mannings from Cross-Section
## Max/Min stage
N2_max_stage, N2_min_stage = PT2['stage(cm)'].max(), PT2['stage(cm)'].min()
## Max stage
#Mannings_Q_from_CrossSection(datadir+'Q/Cross_Section_Surveys/N2_cross_section.xlsx','N2-1_m',Slope=N2_Slope,Manning_n=N2_n,k=N2_k,stage_start=N2_max_stage/100,show=True,save=True,filename=datadir+'Q/Cross_Section_Surveys/N2_Mannings_Cross_Section_stage_max')
## Min stage
#Mannings_Q_from_CrossSection(datadir+'Q/Cross_Section_Surveys/N2_cross_section.xlsx','N2-1_m',Slope=N2_Slope,Manning_n=N2_n,k=N2_k,stage_start=N2_min_stage/100,show=True,save=True,filename=datadir+'Q/Cross_Section_Surveys/N2_Mannings_Cross_Section_stage_min')

## Read N2_Man Discharge from .csv, or calculate new if needed
if 'N2_Man' not in locals():
    try:
        print 'Loading Mannings Q for DAM from CSV'
        N2_Man_reduced = pd.DataFrame.from_csv(datadir+'Q/Manning_Q_files/N2_Man_reduced.csv')
        N2_Man = pd.DataFrame.from_csv(datadir+'Q/Manning_Q_files/N2_Man.csv')
    except:
        print 'Calculate Mannings Q for N2 and saving to CSV'
        N2_stage_reduced = Nuuuli_stage_data['N2'].dropna().round(0).drop_duplicates().order()
        N2_Man_reduced = Mannings_Q_from_stage_data(datadir+'Q/Cross_Section_Surveys/N2_cross_section.xlsx','N2-1_m',Slope=N2_Slope,Manning_n=N2_n,k=N2_k,stage_data=N2_stage_reduced)
        N2_Man_reduced.to_csv(datadir+'Q/Manning_Q_files/N2_Man_reduced.csv')
        N2_stage= Nuuuli_stage_data['N2']
        N2_Man= Mannings_Q_from_stage_data(datadir+'Q/Cross_Section_Surveys/N2_cross_section.xlsx','N2-1_m',Slope=N2_Slope,Manning_n=N2_n,k=N2_k,stage_data=N2_stage)
        N2_Man.to_csv(datadir+'Q/Manning_Q_files/N2_Man.csv')
        pass
#N2_Man_r2 = Manning_AV_r2(N2_Man_reduced,N2stageDischarge)
N2_Man_r2 = 'N/A'

### Compare Discharge Ratings from different methods
def plotQratingN2(ms=6,show=False,log=False,save=False,filename=figdir+''): ## Rating Curves
    mpl.rc('lines',markersize=ms)
    title="Water Discharge Ratings for N2"
    fig, (site_N2, site_N2_zoom)= plt.subplots(1,2,figsize=(8,4))
    xy = np.linspace(0,8000,8000)
    site_N2.text(0.1,0.95,'(a)',verticalalignment='top', horizontalalignment='right',transform=site_N2.transAxes,color='k',fontsize=10,fontweight='bold')
    site_N2_zoom.text(0.1,0.95,'(b)',verticalalignment='top', horizontalalignment='right',transform=site_N2_zoom.transAxes,color='k',fontsize=10,fontweight='bold')
    #N2 AV Measurements and Rating Curve
    site_N2.plot(N2stageDischarge['Q-AV(L/sec)'][start2013:stop2013],N2stageDischarge['stage(cm)'][start2013:stop2013],'^',color='k',fillstyle='none',label='AV 2013') 
    site_N2.plot(N2stageDischarge['Q-AV(L/sec)'][start2014:stop2014],N2stageDischarge['stage(cm)'][start2014:stop2014],'s',color='k',fillstyle='none',label='AV 2014') 
    site_N2_zoom.plot(N2stageDischarge['Q-AV(L/sec)'][start2013:stop2013],N2stageDischarge['stage(cm)'][start2013:stop2013],'^',color='k',fillstyle='none',label='AV 2013') 
    site_N2_zoom.plot(N2stageDischarge['Q-AV(L/sec)'][start2014:stop2014],N2stageDischarge['stage(cm)'][start2014:stop2014],'s',color='k',fillstyle='none',label='AV 2014') 

    ## N2 MODELS
    ## Mannings for AV measurements
    site_N2.plot(N2stageDischarge['Q-AManningV(L/sec)'],N2stageDischarge['stage(cm)'],'.',ls='None',c='grey',label='A-ManningsV')
    site_N2_zoom.plot(N2stageDischarge['Q-AManningV(L/sec)'],N2stageDischarge['stage(cm)'],'.',ls='None',c='grey',label='A-ManningsV')
    ## N2 Power
    N2_AVpower = powerfunction(N2stageDischarge['Q-AV(L/sec)'],N2stageDischarge['stage(cm)'])    
    PowerFit(N2stageDischarge['Q-AV(L/sec)'],N2stageDischarge['stage(cm)'],xy,site_N2,c='grey',ls='-',label='AV power law '+r'$r^2$'+"%.2f"%N2_AVpower['r2'])    
    PowerFit(N2stageDischarge['Q-AV(L/sec)'],N2stageDischarge['stage(cm)'],xy,site_N2_zoom,c='grey',ls='-',label='AV power law '+r'$r^2$'+"%.2f"%N2_AVpower['r2'])        
    ## N2 Mannings from stream survey
    N2_ManQ, N2_Manstage = N2_Man_reduced['Q(m3/s)']*1000, N2_Man_reduced['stage(m)']*100
    site_N2.plot(N2_ManQ,N2_Manstage,'-',markersize=2,c='k',label='Mannings: n='+str(N2_n)+r'$ r^2$'+N2_Man_r2)
    site_N2_zoom.plot(N2_ManQ,N2_Manstage,'-',markersize=2,c='k',label='Mannings')
    ## Label point -click
    labelindex_subplot(site_N2, N2stageDischarge.index,N2stageDischarge['Q-AV(L/sec)'],N2stageDischarge['stage(cm)'])
    labelindex_subplot(site_N2_zoom, N2stageDischarge.index,N2stageDischarge['Q-AV(L/sec)'],N2stageDischarge['stage(cm)'])
    ## Label subplots    
    site_N2.set_ylabel('Stage(cm)'),site_N2.set_xlabel('Q(L/sec)'),site_N2_zoom.set_xlabel('Q(L/sec)')
    ## Format subplots
    site_N2.set_ylim(0,PT2['stage(cm)'].max()+10)#,site_N2.set_xlim(0,N2_AVnonLinear(PT1['stage'].max()+10))
    site_N2_zoom.set_ylim(20,60), site_N2_zoom.set_xlim(0,1600)
    ## Legends
    site_N2.legend(loc='lower right',fancybox=True)  
    ## Figure title
    #plt.suptitle(title,fontsize=16)
    fig.canvas.manager.set_window_title('Figure : '+title) 
    logaxes(log,fig)
    for ax in fig.axes:
        ax.autoscale_view(True,True,True)
    plt.tight_layout(pad=0.1)
    show_plot(show,fig)
    savefig(save,filename)
    return
plotQratingN2(ms=6,show=True,log=False,save=True,filename=figdir+'Stage-Q rating N2')
#plotQratingN2(ms=6,show=True,log=True,save=False,filename=figdir+'')



















