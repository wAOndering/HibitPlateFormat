'''
2022-05-19 
contact: vaissire.t@ufl.edu

This is a set of function to facilitate downstream analysis and vizualization of assays
performed in a 384 format.
The format is flexible however some of the function would need to be adapted for that
especially the get file function

Important all the STD is based on the STD from the 

## TODO: THIS IS A MUST  + CLEAN UP includes some of the method in Class
## TODO: try to set this up with input instead or argparse
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
import copy
import matplotlib
from pathlib import Path
import argparse
import sys
import subprocess
try:
    import pyfiglet
except:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pyfiglet'])
    import pyfiglet
#### LOAD set of custom functions required for further processing
def quickConversion(tmp, myCol=None, option=1):
    if option ==0:
        tmp = tmp.reset_index()
    elif option ==1:
        tmp = tmp.reset_index()
        if tmp.columns.nlevels > 1:
            tmp.columns = ['_'.join(col) for col in tmp.columns] 
        tmp.columns = tmp.columns.str.replace('[_]','')
        if myCol:
            tmp = tmp.rename(columns={np.nan: myCol})
    return tmp

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
    combo['rowCol'] = combo['row'].astype(str)+'_'+combo['col'].astype(str)

    return combo

def getPlateFormat(comboData, myVal, mPath):
    tmp = pd.pivot_table(comboData, values=myVal, index=['row'], columns=['col'])
    tmp.to_csv(mPath+os.sep+'output'+os.sep+'normPlate_'+myVal+'.csv')
    tmp = tmp.to_numpy()

    return tmp

def getListOutofLim(combo, hexColor, mPath):
    '''
    get the list of elements outside the bos 
    '''
    ### geth the value for the red 3SD data 
   
    saveName = mPath+os.sep+'output'+os.sep
    ###---------------------------------------------------------------------------------
    xMark = stdDevbound(combo.loc[combo['sample']=='DMSO','norm_ffluc']) #combo['norm_ffluc'] - original value to look at SD across the entire plate
    yMark = stdDevbound(combo.loc[combo['sample']=='DMSO','norm_nanoluc'])#combo['norm_nanoluc']  - original value to look at SD across the entire plate
    ###---------------------------------------------------------------------------------
    
    ### substraction to get everyting EXCEPT the center areas
    xmask = (combo['norm_ffluc']>xMark[0]) & (combo['norm_ffluc']<xMark[1])
    ymask = (combo['norm_nanoluc']>yMark[0]) & (combo['norm_nanoluc']<yMark[1])
    tmp = combo[~(xmask & ymask)]


    ## assign different area for downstream plotting for samples out of the limit
    tmpCat = copy.copy(combo)
    tmpCat.loc[(combo['norm_nanoluc']<yMark[0]) & (combo['norm_ffluc']>xMark[1]),'zone'] = 2
    tmpCat.loc[(combo['norm_nanoluc']<yMark[0]) & (combo['norm_ffluc']<xMark[1]),'zone'] = 1
    tmpCat.loc[(combo['norm_nanoluc']<yMark[0]) & (combo['norm_ffluc']<xMark[0]),'zone'] = 0

    tmpCat.loc[(combo['norm_nanoluc']>yMark[1]) & (combo['norm_ffluc']>xMark[1]),'zone'] = 8
    tmpCat.loc[(combo['norm_nanoluc']>yMark[1]) & (combo['norm_ffluc']<xMark[1]),'zone'] = 7
    tmpCat.loc[(combo['norm_nanoluc']>yMark[1]) & (combo['norm_ffluc']<xMark[0]),'zone'] = 6

    tmpCat.loc[ymask & (combo['norm_ffluc']>xMark[1]),'zone'] = 5
    tmpCat.loc[ymask & (combo['norm_ffluc']<xMark[1]),'zone'] = 4
    tmpCat.loc[ymask & (combo['norm_ffluc']<xMark[0]),'zone'] = 3


    ### substraction to get everyting EXCEPT the center areas
    figure_mosaic =  '''
    AABC
    AABC
    '''

    fig, axes = plt.subplot_mosaic(mosaic=figure_mosaic, figsize=(11,5))
    #for label, ax in axes.items():
    #   ax.text(0.1,0.1, label, color='blue')

    for k in [1/3,2/3]:
        axes['A'].axvline(k, linestyle='dashed', color='red')
        axes['A'].axhline(k, linestyle='dashed', color='red')
    axes['A'].set_ylabel('nanoluc')
    axes['A'].set_xlabel('ffluc')

    ### enumerate a list 
    tmpall = []
    lst = [1/6, 3/6, 5/6]
    for k in lst:
        for j in lst:
            tmpl = [j,k]
            tmpall.append(tmpl)

    for i, j in enumerate(tmpall):
        # print(i, j)
        axes['A'].text(j[0],j[1], 'zone'+str(i), color='blue')


    census = quickConversion(tmpCat.groupby(['sample','zone']).agg({'zone':['count']}))
    censusTot = quickConversion(census.groupby(['sample']).agg({'zonecount':['sum']}))
    census = pd.merge(census, censusTot, on = ['sample'])
    census['ratio'] = census['zonecount']/census['zonecountsum']

    sns.barplot(x="zone", y="zonecount", hue="sample", data=census, ax = axes['B'], palette=hexColor)
    sns.barplot(x="zone", y="ratio", hue="sample", data=census, ax = axes['C'], palette=hexColor)
    plt.tight_layout()
    plt.savefig(saveName+'samplesCounts.png')
    plt.savefig(saveName+'samplesCounts.pdf')
    plt.close('all')

    return tmp


def getColorListHex(mypalette='Paired', dmsoCol = None, n=12):
        hexColor = [dmsoCol]
        for i in range(n):
            if i ==0:
                continue
            dmsoColtmp = matplotlib.cm.get_cmap(mypalette)(i)
            dmsoColtmp = matplotlib.colors.rgb2hex(dmsoColtmp)
            hexColor.append(dmsoColtmp)
            # print(dmsoColtmp)

        return hexColor

def getThePlot(combo, mPath, style='scatter', allPlot = True, logScale=None, dmsoCol = None):
    '''
    logScale: either None for linear or which ever base desired 
    mPath: correspond to the output path

    '''

    print('Generating plots.....')
    fig, ax = plt.subplots(2,2,figsize=([16, 16]))
    lowerlim =  min(axesParam(combo, False)[0][0], axesParam(combo, False)[1][0])# logScale cannot be 0 and should avoid artificial truncation from lower lim
   

    ##################################
    ## Manipulate the color
    ##################################

    if dmsoCol == None:
        dmsoCol = matplotlib.cm.get_cmap('Paired')(0)
        dmsoCol = matplotlib.colors.rgb2hex(dmsoCol)

    hexColor = getColorListHex(mypalette='Paired', dmsoCol = dmsoCol)
    hexColor = hexColor[:len(np.unique(combo['sample']))]
    ##################################
    ## Definition of the graph limit
    ##################################
    if myCustomLimitOpt == str(1):
        myXaxis = [0, 2]
        myYaxis = [0, 2]
    if myCustomLimitOpt == str(2):
        myXaxis = axesParam(combo)[0]
        myYaxis = axesParam(combo)[1]
    if myCustomLimitOpt == str(3):
        myXaxis = float(myCustomLimitOptX_low), float(myCustomLimitOptX_high)
        myYaxis = [float(myCustomLimitOptY_low), float(myCustomLimitOptY_high)]

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
        sns.scatterplot(data=combo, x="norm_ffluc", y="norm_nanoluc", hue="sample", alpha=0.9, palette=hexColor, ax=ax[1][0])
        sns.scatterplot(data=combo[combo['sample']=='DMSO'], x="norm_ffluc", y="norm_nanoluc", alpha=0.9, color=dmsoCol,ax=ax[1][1])

    elif style == 'kde':
        sns.kdeplot(data=combo, x="norm_ffluc", y="norm_nanoluc", hue="sample", alpha=0.9, palette=hexColor, ax=ax[1][0])
        sns.kdeplot(data=combo[combo['sample']=='DMSO'], x="norm_ffluc", y="norm_nanoluc", alpha=0.9, color=dmsoCol,ax=ax[1][1])

    ### geth the value for the red 3SD data 
    ###---------------------------------------------------------------------------------
    xMark = stdDevbound(combo.loc[combo['sample']=='DMSO','norm_ffluc']) #combo['norm_ffluc'] - original value to look at SD across the entire plate
    yMark = stdDevbound(combo.loc[combo['sample']=='DMSO','norm_nanoluc'])#combo['norm_nanoluc']  - original value to look at SD across the entire plate
    ###---------------------------------------------------------------------------------

    ### introduction of log scale
    ###---------------------------------------------------------------------------------
    if logScale != None:
        ax[1][0].set_xscale("log", base=logScale)
        ax[1][0].set_yscale("log", base=logScale)
        ax[1][1].set_xscale("log", base=logScale)
        ax[1][1].set_yscale("log", base=logScale)
        saveName = mPath+os.sep+'output'+os.sep+'figSummary_log'+str(logScale)+'_'+style
        
        myXaxis = [min(0.2, lowerlim), myXaxis[1]]
        myYaxis = [min(0.2, lowerlim), myYaxis[1]]   
    else:
        saveName = mPath+os.sep+'output'+os.sep+'figSummary_'+style        
    ###---------------------------------------------------------------------------------
    
    for i in ax[1]:
        i.set_xlim(myXaxis)
        i.set_ylim(myYaxis)
        for j in [0,1]:
            # print(xMark[j])
            i.axvline(xMark[j], linestyle='dashed', color='red')
            i.axhline(yMark[j], linestyle='dashed', color='red')

    ax[0][0].set_title('Plate map of ffluc', fontsize= 12)
    ax[0][1].set_title('Plate map of nanoluc', fontsize= 12)
    ax[1][0].set_title('Scatter for all the samples (red lines: 3*SD of DMSO)', fontsize= 12)
    ax[1][1].set_title('Scatter for DMSO only (red lines: 3*SD of DMSO)', fontsize= 12)

    plt.savefig(saveName+'.png')
    plt.savefig(saveName+'.pdf')
    plt.close('all')
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

            for k in [0,1]:
                # print(xMark[j])
                i.axvline(xMark[k], linestyle='dashed', color='red')
                i.axhline(yMark[k], linestyle='dashed', color='red')
                # i.set_title('Sample: '+str(j), fontsize= 12)
                i.text(0.1, 0.2, 'Sample: '+str(j), fontsize= 12)

                ### introduction of log scale
                ###---------------------------------------------------------------------------------
                if logScale != None:
                    i.set_xscale("log", base=logScale)
                    i.set_yscale("log", base=logScale)
                    myXaxis = [min(0.2, lowerlim), myXaxis[1]]
                    myYaxis = [min(0.2, lowerlim), myYaxis[1]]  
                    saveName = mPath+os.sep+'output'+os.sep+'figSummaryInd_log'+str(logScale)
                else:
                    saveName = mPath+os.sep+'output'+os.sep+'figSummaryInd_'
                ###---------------------------------------------------------------------------------

            i.set_xlim(myXaxis)
            i.set_ylim(myYaxis)

        plt.savefig(saveName+'.png')
        plt.savefig(saveName+'.pdf')
        plt.close('all')
        
        ##################################
        ## Get the outlier samples
        ##################################
        fig, ax = plt.subplots(figsize=([16, 16]))

        for k in [0,1]:
            # print(xMark[j])
            ax.axvline(xMark[k], linestyle='dashed', color='red')
            ax.axhline(yMark[k], linestyle='dashed', color='red')
        
        tmp = getListOutofLim(combo, hexColor, mPath)
        sns.scatterplot(data=combo, x="norm_ffluc", y="norm_nanoluc", s=80, color='black', alpha=0.2)
        sns.scatterplot(data=tmp, x="norm_ffluc", y="norm_nanoluc", s=80, color='red', alpha=0.8)
        for i, j in tmp.iterrows():
            mytext = str(j['sample'])+': '+j['row']+str(j['col'])
            ax.text(j['norm_ffluc'], j['norm_nanoluc'], mytext, fontsize = 12)

        print('')
        print('Table for data falling out of the 3*SD limits - 3SD from the DMSO group')
        print(tmp)

        ### introduction of log scale
        ###---------------------------------------------------------------------------------
        if logScale != None:
            ax.set_xscale("log", base=logScale)
            ax.set_yscale("log", base=logScale)
            myXaxis = [min(0.2, lowerlim), myXaxis[1]]
            myYaxis = [min(0.2, lowerlim), myYaxis[1]]  
            saveName = mPath+os.sep+'output'+os.sep+'figSummaryOuterSample_log'+str(logScale)
        else:
            saveName = mPath+os.sep+'output'+os.sep+'figSummaryOuterSample_'
        ###---------------------------------------------------------------------------------
        ax.set_xlim(myXaxis)
        ax.set_ylim(myYaxis)

        tmp.to_csv(mPath+os.sep+'output'+os.sep+'outerSample.csv')
        
        plt.savefig(saveName+'.png')
        plt.savefig(saveName+'.pdf')
        plt.close('all')

def getCV(mPath, combo, dmsoCol = None):
    '''
    Function to generate and display CV by group 
    CV being defined as std/mean
    '''
    # def custCV(x):
    #     return np.std(x)/np.mean(x)
    saveName = mPath+os.sep+'output'+os.sep

    meltCombo = combo.melt(['row','col','sample', 'rowCol', 'plate'])
    tmp = meltCombo.groupby(['sample','variable']).agg({'value':[min, max, sum, 'count', np.mean, np.std]})
    a = quickConversion(tmp)
    a['CV'] = a['valuestd']/a['valuemean']
    a.to_csv(saveName+'Samples_stats.csv')

    ##################################
    ## Manipulate the color
    ##################################

    if dmsoCol == None:
        dmsoCol = matplotlib.cm.get_cmap('Paired')(0)
        dmsoCol = matplotlib.colors.rgb2hex(dmsoCol)

    hexColor = getColorListHex(mypalette='Paired', dmsoCol = dmsoCol)
    hexColor = hexColor[:len(np.unique(combo['sample']))]
    fig, ax = plt.subplots(figsize=(3,6))

    tmpCv = a[a['variable'].str.contains('norm')] # plot CV only for the normalized data 
    sns.barplot(x="variable", y="CV", hue="sample", data=tmpCv, ax = ax, palette=hexColor)
    ax.set_ylim(0,0.5)
    plt.tight_layout()
    plt.savefig(saveName+'CVplot.png')
    plt.savefig(saveName+'CVplot.pdf')
    plt.close('all')

## default or set user inputs criteria
outlierCut = 1

## get the arguments after the file 
# mPath = sys.argv[1]
# if os.path.exists(mPath):
#     print ("Folder exist")
#     os.makedirs(mPath+os.sep+'output',exist_ok=True)

# fileType = input('Enter the file type to be working with options (csv or xlsx):')
print('#########################################################################################################')
print(pyfiglet.figlet_format("HiBit",font='isometric3', width=10000))
print('#########################################################################################################')



print('')
print('-------------------------------------------')

addPlate = '1'
folderList = []
while addPlate == str(1):
    print('Select the FOLDER containing the files to analyze')

    tmpFol = input("Drag the FOLDER and press Enter:")
    tmpFol = tmpFol.replace('\\','/')
    tmpFol = tmpFol.replace('"','')
    folderList.append(tmpFol)
    print('')
    print('Do additional plates need to be included and averaged with the previous one(s)')
    addPlate = input("Answer (1: yes or 0: no): ")

print('-------------------------------------------')
print('')


print('')
print('-------------------------------------------')
print('Select the color for the DMSO control group:')
print("1: default orange ('#EC9540')")
print("2: Brighter orange ('#F05A28')")
print("3: enter the HEX value of the desired color Brighter orange ('#F05A28')")
dmsoCol = input('make a selection (1, 2 or #F05A28):')
if dmsoCol == str(1):
    dmsoCol = '#EC9540'
elif dmsoCol == str(2):
    dmsoCol = '#F05A28'

print('')
print("FYI:" )
print("color options are defaulted to ['Paired'] in color maps (link 1 below) for more color options visit:" )
# print("in the future create a custom color map generation where user select colors then become of cmap and can be substituded from Paired")
print('link 1: https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html')
print('link 2: https://colorbrewer2.org/')
print('link 3: https://color.adobe.com/')

# input('Press Enter to continue')
print('-------------------------------------------')
print('')


print('')
print('-------------------------------------------')
print('Select output types:')
print('1: option1 - 2 groups // with 1 group DMSO on the first and last 2 columns vs 1 sample in between those wells')
print('2: option1 - 12 groups // with 1 group DMSO on the first and last 2 columns and other samples being duplicated columns')
customSampleFormat = input('make a selection (1 or 2):')
print(customSampleFormat)
print('-------------------------------------------')
print('')

print('-------------------------------------------')
print('Select the limits to be displayed on the graph options:')
print('1: (default) option with limit set to 2 for both x and y axes')
print("2: have the axes be fit to the data from 0 to 10'%' of maximum range for both x and y axes" )
print("3: input your own axes value for x and y" )
myCustomLimitOpt = input('Make option selection based on what is described above (eg. 1):')
if myCustomLimitOpt == str(3):
    myCustomLimitOptX_low = input('Input the MIN value of x (eg: 0):')
    myCustomLimitOptX_high = input('Input the MAX value of x (eg: 0):')
    myCustomLimitOptY_low = input('Input the MIN value of y (eg: 0):')
    myCustomLimitOptY_high = input('Input the MAX value of y (eg: 0):')
print('-------------------------------------------')
print('')

## this section can be turned back on to have optional graphing

# print('-------------------------------------------')
# print('Graph will be generated on linear scale and also log scale')
# mylogScale = input('Select the base for the log scale (eg. 2 or 10):')
# mylogScale = int(mylogScale)
# print('-------------------------------------------')
# print('')

# mPath = input('Drag the folder containing the file to analyze:')
# print(os.path.exists(mPath))
print(folderList)



#####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### the strategy here is to 
'''
1. generate a full combined plate 
2. an average plate 
3. tracitional separate plates 
'''
#####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# print(glob.glob(mPath+os.sep+'*.*'))
combPlate = []
for i, j in enumerate(folderList):
    os.makedirs(j+os.sep+'output',exist_ok=True)

    print(i,j)
    combo = combineData(j)
    combo['plate'] = i

    combo.to_csv(j+os.sep+'output'+os.sep+'combinedOutput.csv')
    getThePlot(combo, mPath=j, style='scatter', dmsoCol=dmsoCol)
    getThePlot(combo, mPath=j, style='kde', allPlot = False, dmsoCol=dmsoCol)

    # get plot on log 2
    mylogScale = 2
    getThePlot(combo, mPath=j, style='scatter', logScale=mylogScale, dmsoCol=dmsoCol)
    getThePlot(combo, mPath=j, style='kde', allPlot = False, logScale=mylogScale, dmsoCol=dmsoCol)

    # get plot on log 10
    # mylogScale = 10
    # getThePlot(mPath, 'scatter', logScale=mylogScale)
    # getThePlot(mPath, style='kde', allPlot = False, logScale=mylogScale)
    getCV(j, combo)

    combPlate.append(combo)


#####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### the strategy here is to 
'''
1. generate a full combined plate 
2. an average plate 
3. tracitional separate plates 
'''
#####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if len(folderList)>1:
    # this plate can be used to get the combined out put
    combPlate = pd.concat(combPlate)
    # combPlate['locTag'] = combPlate['rowCol']+'_'+combPlate['plate'].astype(str)

    # this plate can be used as the average plate
    combPlateM = pd.melt(combPlate, id_vars = ['row', 'col', 'sample', 'rowCol', 'plate'], 
        value_vars = ['value_ffluc', 'norm_ffluc', 'value_nanoluc', 'norm_nanoluc'])
    combPlateM_tmp = combPlateM.groupby(['row', 'col', 'sample', 'rowCol', 'variable']).agg({'value':np.mean})
    combPlateM_tmp = quickConversion(combPlateM_tmp)
    combPlateM_tmp = pd.pivot_table(combPlateM_tmp, index=['row', 'col', 'sample', 'rowCol'], values='value', columns=['variable'])
    combPlateM_tmp = quickConversion(combPlateM_tmp, option=0)


    ##################################################
    ##### aggregated output 
    ##################################################
    # for ploting this specific case this plate can be used to get the combined out put
    aggOut = folderList[0]+os.sep+'output_aggregated' # for aggregated output file path
    os.makedirs(aggOut,exist_ok=True)
    os.makedirs(aggOut+os.sep+'output',exist_ok=True)
    combPlate['plate'] = 0 # remake the palte to 0
    combPlate = combPlate.reset_index(drop=True)
    combPlate.to_csv(aggOut+os.sep+'combinedOutput.csv')
    getThePlot(combPlate, mPath=aggOut, style='scatter', dmsoCol=dmsoCol)
    getThePlot(combPlate, mPath=aggOut, style='kde', allPlot = False, dmsoCol=dmsoCol)

    # get plot on log 2
    mylogScale = 2
    getThePlot(combPlate, mPath=aggOut, style='scatter', logScale=mylogScale, dmsoCol=dmsoCol)
    getThePlot(combPlate, mPath=aggOut, style='kde', allPlot = False, logScale=mylogScale, dmsoCol=dmsoCol)

    # get plot on log 10
    # mylogScale = 10
    # getThePlot(mPath, 'scatter', logScale=mylogScale)
    # getThePlot(mPath, style='kde', allPlot = False, logScale=mylogScale)
    getCV(aggOut, combPlate)


    ##################################################
    ##### averaged output 
    ##################################################
    avgOut = folderList[0]+os.sep+'ouput_averaged'
    os.makedirs(avgOut,exist_ok=True)
    os.makedirs(avgOut+os.sep+'output',exist_ok=True)
    combPlateM_tmp.to_csv(avgOut+os.sep+'combinedOutput.csv')
    getThePlot(combPlateM_tmp, mPath=avgOut, style='scatter', dmsoCol=dmsoCol)
    getThePlot(combPlateM_tmp, mPath=avgOut, style='kde', allPlot = False, dmsoCol=dmsoCol)

    # get plot on log 2
    mylogScale = 2
    getThePlot(combPlateM_tmp, mPath=avgOut, style='scatter', logScale=mylogScale, dmsoCol=dmsoCol)
    getThePlot(combPlateM_tmp, mPath=avgOut, style='kde', allPlot = False, logScale=mylogScale, dmsoCol=dmsoCol)

    # get plot on log 10
    # mylogScale = 10
    # getThePlot(mPath, 'scatter', logScale=mylogScale)
    # getThePlot(mPath, style='kde', allPlot = False, logScale=mylogScale)
    combPlateM_tmp['plate'] = 0
    getCV(avgOut, combPlateM_tmp)




print('')
print('')
print('The outputs for the individual plates are located in their respective folders')
for i in folderList:
    print('     *'+i+os.sep+'output')

if len(folderList)<1:
    print('The ouputs for the averaged or aggregated plate are save here')
    print('     *'+folderList[0]+os.sep+'ouput_averaged')
    print('     *'+folderList[0]+os.sep+'output_aggregated')

print('')
print('#########################################################################################################')
print(pyfiglet.figlet_format("Bye Bye!",font='isometric3', width=10000))
print('#########################################################################################################')
print('')

input('press Sapce and Enter to exit')


# 'C:/Users/Windows/Desktop/CamiloTest'