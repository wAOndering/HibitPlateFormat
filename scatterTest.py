'''
2022-05-19 
contact: vaissire.t@ufl.edu

This is a set of function to facilitate downstream analysis and vizualization of assays
performed in a 384 format.
The format is flexible however some of the function would need to be adapted for that
especially the get file function

Important all the STD is based on the STD from the 

## TODO: try to set this up with input instead or argparse
## TODO: includes some of the method in Class
## TODO: integrate some options
## TODO: comments the code better
## TODO: update github
## TODO: create a custom color map generation where user select colors then become of cmap and can be substituded from Paired
## TODO: can create a more interactive interface with the command prompt
## TODO: substitution of sample generate by user input samples
## TODO: sample color can be used and determine to create a color map of the plate format
'''


#### LOAD set of python libraries required for further processing
import warnings
warnings.simplefilter("ignore")
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import string
import matplotlib
from pathlib import Path
import argparse
import sys


#### LOAD set of custom functions required for further processing
def createAplate():
    '''
    Function to enable the creation of a plate format 
    '''
    w_dim = 24
    h_dim = 16
    # make an empty data set
    data = np.zeros((h_dim, w_dim))

    # get the control data plate
    data[2:14,0:2] = 1
    data[2:14,22:24] = 1
    
    # get the range for sample
    data[:,2:22]=2

    # create a colormap
    my_cmap = matplotlib.colors.ListedColormap(['#f7f7f7','#ef8a62','#67a9cf'])

    # generate the figure
    fig, ax = plt.subplots(figsize=[10,4])
    mymap = ax.pcolormesh(data, cmap=my_cmap, edgecolors='grey', linewidth=0.5)
    
    # This function formatter will replace integers with target names
    formatter = plt.FuncFormatter(lambda val, loc: ['NA', 'DMSO ctl', 'Samples'][val])
    plt.colorbar(mymap, ticks=[0, 1, 2], format=formatter)
    mymap.set_clim(-0.5, 2.5)     # Set the clim so that labels are centered on each block

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False) #0.5 enable to set the axis in the middle
    ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False) #0.5 enable to set the axis in the middle
    ax.invert_yaxis()

    #labels
    column_labels = [str(x) for x in list(np.arange(24)+1)]
    row_labels = list(string.ascii_uppercase)[:h_dim] 
    ax.set_xticklabels(column_labels, minor=False)
    ax.set_yticklabels(row_labels, minor=False)
    plt.title('Standard plate format')

    plt.savefig('plateFormat.png')
    
def getFile(fileName, dattype):
    '''
    Main function initially set up as default by reading the xlsx file
    which based on pandas version can be a source of error. This has been implemented 
    
    Args:
        fileName (str): name of the file of interest
        dattype (str): name of data to work on either hibit of ffluc, only serve to label column properly
        which is important for downstream merge operations 


    to work on the csv file format as well.
    1- format as long format
    2- get the baseline signal based on DMSO well correpsonding to the frist 2 and last 2 columns of the plate
    minus the 4 wells 
    '''

    if fileName.split(os.sep)[-1].split('.')[-1] == 'xlsx':
        ## reformat and reshape the data to a long format
        tmp = pd.read_excel(fileName, sheet_name='Results', usecols='F:AC', skiprows=9) # read the excel file based on location on the actual data
        plateH = np.shape(tmp)[0] # get the format of the plate height and 
        plateW = np.shape(tmp)[1] # get the format of the plate width
        tmp['row'] = list(string.ascii_uppercase)[:plateH] # create row labels
        tmpM = tmp.melt('row')
        tmpM = tmpM.drop(columns = ['variable']) # remove the variable columns which will be replaced by better formation of the col
        tmpM['col'] = np.repeat(range(plateW),plateH)+1 # create col labels in python this is +1 as 0 is the start
    
    elif fileName.split(os.sep)[-1].split('.')[-1] == 'csv':
        tmpO = pd.read_csv(fileName) #temporary store the csv files prior to resoncstruction
        tmpM = tmpO['WellPosition'].str.split(':', expand=True)
        tmpM.columns = ['row', 'col']
        tmpM['value'] = tmpO['RLU']
        tmpM['col'] = tmpM['col'].astype(int) # the column type when csv is imported is an object need comversion for porper sorting 
    
    # the way the formating is done is when the plate is sorted first by row then by columns'
    # this is defautl for xlsx but adapted for csv
    tmpM = tmpM.sort_values(['col','row'])
    tmpM = tmpM.reset_index(drop=True)


    ##------------------------------------------------------------------
    ## choose the type of labeling customSampleFormat
    ##------------------------------------------------------------------

    # this line is important as it determines the plate format and assume that the samples are repeated in duplicated columns
    tmpM['sample'] = np.repeat(range(int(plateW/2)),plateH*2)+1 
    # DMSO control is established as the set of 2 columns
    tmpM.loc[tmpM['sample'].isin([1,12]),'sample'] ='DMSO'
    tmpM = tmpM[['row','col','sample','value']]
    print(customSampleFormat)
    if customSampleFormat == '1':
        tmpM.loc[~(tmpM['sample']=='DMSO'),'sample'] = 'EXP'

    # get all the DMSO sample defined in the sample section see above
    # and take out the wells which are A, B, O, P
    dmsoCtl = tmpM[(tmpM['sample'].isin(['DMSO'])) & (~tmpM['row'].isin(['A','B','O','P']))]

    
    lwr_b = stdDevbound(dmsoCtl['value'])[0]
    hir_b = stdDevbound(dmsoCtl['value'])[1]
    dmsoCtlwoOut = dmsoCtl[(dmsoCtl['value']>lwr_b) & (dmsoCtl['value']<hir_b)]
    outlier_def = dmsoCtl[~(dmsoCtl['value']>lwr_b) & (dmsoCtl['value']<hir_b)]
    
    print('For the file:')
    # print(fileName)
    print('The outlier range is: ', [int(lwr_b),int(hir_b)])
    print('')
    print('Outliers within the DMSO are listed below:')
    print(outlier_def)
    print('')
    print('Outliers are removed by default option can be set to keep them or change their definition')
    print('Outliers will be displayed as white box in the heatmap')
    if outlierCut == 1:
        tmpMnoOut = tmpM.drop(outlier_def.index, axis = 0)
    print('')
    print('------------------------------------------------------')

    tmpMnoOut['norm'] = tmpMnoOut['value']/np.mean(dmsoCtlwoOut['value'])
    tmpMnoOut = tmpMnoOut.rename(columns={'value':'value_'+dattype,'norm':'norm_'+dattype})

    return tmpMnoOut

def axesParam(data, atZero=True):
    '''
    This function returns the limit for the graph. When atZero: default(True) 
    this will set y and x axis at zero with common maxima based on the maximal value
    in both the norm_nanoluc or norm_ffluc  data sets

    data: correpsond to the combine data file generate combo
    atZero: default(True) 
    '''
    
    myXlim = [min(data['norm_ffluc'])*0.9, max(data['norm_ffluc'])*1.1]
    myYlim = [min(data['norm_nanoluc'])*0.9, max(data['norm_nanoluc'])*1.1]

    if atZero == True:
        myXlim = [0,1.1*max(max(data['norm_nanoluc']),max(data['norm_ffluc']))]
        myYlim = [0,1.1*max(max(data['norm_nanoluc']),max(data['norm_ffluc']))]

    return myXlim, myYlim

def stdDevbound(dataArray, custSD=3):
    '''
    Enables to obtain lower and higher bounds of the data based on the criteria of 3 standard deviation
    ##TODO within the class when created append the thresholding methods to the filename

    dataArray: correspond to the array for which we want to obtain the limit of interest
    '''
    outlierCutoffValue = custSD*np.std(dataArray)
    lwr_b = np.mean(dataArray)-outlierCutoffValue
    hir_b = np.mean(dataArray)+outlierCutoffValue

    return [lwr_b, hir_b]

def combineData(mPath, fileType = 'xlsx'):
    '''
    Function which enable to combine the ffluc and nanoluc data within the same file
    the initial fileType developped on was xlsx
    '''
    dattype = 'ffluc'
    # here the set are present to be able to deal with potential presence of temporary open excel files ~ which can create problem
    fflname = getFile(list(set(glob.glob(mPath+os.sep+'*'+dattype+'*.'+fileType))-set(glob.glob(mPath+os.sep+'~*'+dattype+'*.'+fileType)))[0], dattype = dattype)

    dattype = 'nanoluc'
    hibname = getFile(list(set(glob.glob(mPath+os.sep+'*'+dattype+'*.'+fileType))-set(glob.glob(mPath+os.sep+'~*'+dattype+'*.'+fileType)))[0], dattype = dattype)

    combo = pd.merge(fflname,hibname,on=['row','col','sample'])

    return combo

def getPlateFormat(comboData, myVal, mPath):
    tmp = pd.pivot_table(comboData, values=myVal, index=['row'], columns=['col'])
    tmp.to_csv(mPath+os.sep+'output'+os.sep+'normPlate_'+myVal+'.csv')
    tmp = tmp.to_numpy()

    return tmp

def getListOutofLim(combo):
    '''
    get the list of elements outside the bos 
    '''
    ### geth the value for the red 3SD data 
    ###---------------------------------------------------------------------------------
    xMark = stdDevbound(combo.loc[combo['sample']=='DMSO','norm_ffluc']) #combo['norm_ffluc'] - original value to look at SD across the entire plate
    yMark = stdDevbound(combo.loc[combo['sample']=='DMSO','norm_nanoluc'])#combo['norm_nanoluc']  - original value to look at SD across the entire plate
    ###---------------------------------------------------------------------------------
    tmp = combo[~((combo['norm_ffluc']>xMark[0]) & (combo['norm_ffluc']<xMark[1]))]


    xmask = (combo['norm_ffluc']>xMark[0]) & (combo['norm_ffluc']<xMark[1])
    ymask = (combo['norm_nanoluc']>yMark[0]) & (combo['norm_nanoluc']<yMark[1])
    tmp = combo[~(xmask & ymask)]

    return tmp



    tmp = tmp[~((tmp['norm_nanoluc']>yMark[0]) & (tmp['norm_ffluc']<yMark[1]))]
    combo[combo['norm_ffluc']]

def getThePlot(mPath, style='scatter', allPlot = True):
    print('Generating plots.....')

    fig, ax = plt.subplots(2,2,figsize=([16, 16]))

    combo = combineData(mPath)
    ##################################
    ## Get the plate map HEATMAP
    ##################################
    ## generate heatmaps for both the ffluc and nanoluc

    myVmin = min(np.array(axesParam(combo, False)).flatten())
    myVmax = max(np.array(axesParam(combo, False)).flatten())
    ax[0][0].imshow(getPlateFormat(combo, 'norm_ffluc', mPath), interpolation = 'nearest', aspect='auto', vmin=myVmin, vmax=myVmax) # cmap= 'bwr')
    ax[0][1].imshow(getPlateFormat(combo, 'norm_nanoluc', mPath), interpolation = 'nearest', aspect='auto', vmin=myVmin, vmax=myVmax) # cmap= 'bwr')
    
    for i in ax[0]:
        i.set_xlabel('plate columns')
        i.set_ylabel('plate rows')

    ##################################
    ## Get the scatter
    ##################################
    if style == 'scatter':
        sns.scatterplot(data=combo, x="norm_ffluc", y="norm_nanoluc", hue="sample", alpha=0.9, palette='Paired', ax=ax[1][0])

        dmsoCol = matplotlib.cm.get_cmap('Paired')(0)
        dmsoCol = matplotlib.colors.rgb2hex(dmsoCol)

        sns.scatterplot(data=combo[combo['sample']=='DMSO'], x="norm_ffluc", y="norm_nanoluc", alpha=0.9, color=dmsoCol,ax=ax[1][1])

    elif style == 'kde':
        sns.kdeplot(data=combo, x="norm_ffluc", y="norm_nanoluc", hue="sample", alpha=0.9, palette='Paired', ax=ax[1][0])

        dmsoCol = matplotlib.cm.get_cmap('Paired')(0)
        dmsoCol = matplotlib.colors.rgb2hex(dmsoCol)

        sns.kdeplot(data=combo[combo['sample']=='DMSO'], x="norm_ffluc", y="norm_nanoluc", alpha=0.9, color=dmsoCol,ax=ax[1][1])

    ### geth the value for the red 3SD data 
    ###---------------------------------------------------------------------------------
    xMark = stdDevbound(combo.loc[combo['sample']=='DMSO','norm_ffluc']) #combo['norm_ffluc'] - original value to look at SD across the entire plate
    yMark = stdDevbound(combo.loc[combo['sample']=='DMSO','norm_nanoluc'])#combo['norm_nanoluc']  - original value to look at SD across the entire plate
    ###---------------------------------------------------------------------------------

    for i in ax[1]:
        i.set_xlim(axesParam(combo)[0])
        i.set_ylim(axesParam(combo)[1])
        for j in [0,1]:
            # print(xMark[j])
            i.axvline(xMark[j], linestyle='dashed', color='red')
            i.axhline(yMark[j], linestyle='dashed', color='red')

    ax[0][0].set_title('Plate map of ffluc', fontsize= 12)
    ax[0][1].set_title('Plate map of nanoluc', fontsize= 12)
    ax[1][0].set_title('Scatter for all the samples (red lines: 3*SD of DMSO)', fontsize= 12)
    ax[1][1].set_title('Scatter for DMSO only (red lines: 3*SD of DMSO)', fontsize= 12)

    saveName = mPath+os.sep+'output'+os.sep+'figSummary_'+style
    plt.savefig(saveName+'.png')
    plt.savefig(saveName+'.pdf')

    ##################################
    ## Get the plot for individual samples
    ##################################
    
    if allPlot == True:
        fig, ax = plt.subplots(4,3,figsize=([16, 16]))

        for idx, (i, j) in enumerate(zip(ax.flatten(), combo['sample'].unique())):
            # print(idx, i, j)
            dmsoCol = matplotlib.cm.get_cmap('Paired')(idx)
            dmsoCol = matplotlib.colors.rgb2hex(dmsoCol)
            sns.scatterplot(data=combo[combo['sample']==j], x="norm_ffluc", y="norm_nanoluc", alpha=0.9, color=dmsoCol,ax=i)
            sns.kdeplot(data=combo[combo['sample']==j], x="norm_ffluc", y="norm_nanoluc", alpha=0.3, color=dmsoCol,ax=i)
            i.set_xlim(axesParam(combo)[0])
            i.set_ylim(axesParam(combo)[1])
            for k in [0,1]:
                # print(xMark[j])
                i.axvline(xMark[k], linestyle='dashed', color='red')
                i.axhline(yMark[k], linestyle='dashed', color='red')
                # i.set_title('Sample: '+str(j), fontsize= 12)
                i.text(0.1, 0.2, 'Sample: '+str(j), fontsize= 12)

        saveName = mPath+os.sep+'output'+os.sep+'figSummaryInd_'
        plt.savefig(saveName+'.png')
        plt.savefig(saveName+'.pdf')

        ##################################
        ## Get the outlier samples
        ##################################
        fig, ax = plt.subplots(figsize=([16, 16]))
        ax.set_xlim(axesParam(combo)[0])
        ax.set_ylim(axesParam(combo)[1])
        for k in [0,1]:
            # print(xMark[j])
            ax.axvline(xMark[k], linestyle='dashed', color='red')
            ax.axhline(yMark[k], linestyle='dashed', color='red')
        
        tmp = getListOutofLim(combo)
        sns.scatterplot(data=combo, x="norm_ffluc", y="norm_nanoluc", s=80, color='black', alpha=0.2)
        sns.scatterplot(data=tmp, x="norm_ffluc", y="norm_nanoluc", s=80, color='red', alpha=0.8)
        for i, j in tmp.iterrows():
            mytext = str(j['sample'])+': '+j['row']+str(j['col'])
            ax.text(j['norm_ffluc'], j['norm_nanoluc'], mytext, fontsize = 12)

        print('')
        print('Table for data falling out of the 3*SD limits - 3SD from the DMSO group')
        print(tmp)

        tmp.to_csv(mPath+os.sep+'output'+os.sep+'outerSample.csv')
        saveName = mPath+os.sep+'output'+os.sep+'figSummaryOuterSample_'
        plt.savefig(saveName+'.png')
        plt.savefig(saveName+'.pdf')

def quickConversion(tmp, myCol=None):
    '''
    tools to quickly convert the output of groupby
    '''
    tmp = tmp.reset_index()
    if tmp.columns.nlevels > 1:
        tmp.columns = ['_'.join(col) for col in tmp.columns] 
    tmp.columns = tmp.columns.str.replace('[_]','')
    if myCol:
        tmp = tmp.rename(columns={np.nan: myCol})
    return tmp

def getCV(mPath, combo):
    '''
    Function to generate and display CV by group 
    CV being defined as std/mean
    '''
    # def custCV(x):
    #     return np.std(x)/np.mean(x)

    meltCombo = combo.melt(['row','col','sample'])
    tmp = meltCombo.groupby(['sample','variable']).agg({'value':[min, max, sum, 'count', np.mean, np.std]})
    a = quickConversion(tmp)
    a['CV'] = a['valuestd']/a['valuemean']
    a.to_csv(mPath+os.sep+'output'+os.sep+'Samples_stats.csv')

## default or set user inputs criteria
outlierCut = 1

## get the arguments after the file 
mPath = sys.argv[1]
if os.path.exists(mPath):
    print ("Folder exist")
    os.makedirs(mPath+os.sep+'output',exist_ok=True)

# fileType = input('Enter the file type to be working with options (csv or xlsx):')
print('')
print('-------------------------------------------')
print('Select output types:')
print('1: option1 - 2 groups // with 1 group DMSO on the first and last 2 columns vs 1 sample in between those wells')
print('2: option1 - 12 groups // with 1 group DMSO on the first and last 2 columns and other samples being duplicated columns')
customSampleFormat = input('make a selection (1 or 2):')
print(customSampleFormat)
# print("color options are defaulted to ['Paired'] in color maps (link 1 below) for more color options visit:" )
# print("in the future create a custom color map generation where user select colors then become of cmap and can be substituded from Paired")
# print('link 1: https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html')
# print('link 2: https://colorbrewer2.org/')
# input('Press Enter to continue')
print('-------------------------------------------')
print('')
# mPath = input('Drag the folder containing the file to analyze:')
# print(os.path.exists(mPath))
print(mPath)
# print(glob.glob(mPath+os.sep+'*.*'))
combo = combineData(mPath)
combo.to_csv(mPath+os.sep+'output'+os.sep+'combinedOutput.csv')
getThePlot(mPath, 'scatter')
getThePlot(mPath, style='kde', allPlot = False)
getCV(mPath, combo)
# 'C:/Users/Windows/Desktop/CamiloTest'