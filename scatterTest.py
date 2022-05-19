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

## TODO: try to set this up with input instead or argparse
## TODO: implement csv vs xlsx options 
## TODO: includes some of the method in Class
## TODO: integrate some options
## TODO: comments the code better
## TODO: update github
## TODO: can create a more interactive interface with the command prompt

def normalizePlate(dat):
    '''DEPRECATED DO NOT USE DEPRECATED DO NOT USE DEPRECATED DO NOT USE 
    DEPRECATED DO NOT USE DEPRECATED DO NOT USE DEPRECATED DO NOT USE 
    '''

    # get the series corresponding to first and last 2 columns 
    avgSeries = dat.iloc[2:14, np.r_[0:2,22,23]].values.flatten()
    
    # takeaway outlier based on the 3 standard deviation
    avgSeries = avgSeries[np.logical_and((avgSeries>np.mean(avgSeries)-3*np.std(avgSeries)), (avgSeries<np.mean(avgSeries)+3*np.std(avgSeries)))]

    avgDMSO = np.mean(avgSeries)

    datNorm = dat/avgDMSO

    colName = np.repeat(range(12),2)
    colName = ['col'+str(x) for x in colName]

    datNorm.columns = colName
    meltDat = datNorm.melt()

    # below 1 correspond to the average as it has been normalized above
    avgSeriesNorm = datNorm.iloc[2:14, np.r_[0:2,22,23]].values.flatten()
    
    # takeaway outlier based on the 3 standard deviation
    avgSeriesNorm = avgSeriesNorm[np.logical_and((avgSeriesNorm>np.mean(avgSeriesNorm)-3*np.std(avgSeriesNorm)), (avgSeriesNorm<np.mean(avgSeriesNorm)+3*np.std(avgSeriesNorm)))]
    stdLow = 1-3*np.std(avgSeriesNorm)
    stdHigh = 1+3*np.std(avgSeriesNorm)

    return meltDat, [stdLow, stdHigh]

def getFile(fileName, dattype):

    if fileName.split(os.sep)[-1].split('.')[-1] == 'xlsx':
        tmp = pd.read_excel(fileName, sheet_name='Results', usecols='F:AC', skiprows=9)
        plateH = np.shape(tmp)[0]
        plateW = np.shape(tmp)[1]
        tmp['row'] = list(string.ascii_uppercase)[:plateH]
        tmpM = tmp.melt('row')
        tmpM = tmpM.drop(columns = ['variable'])
        tmpM['col'] = np.repeat(range(plateW),plateH)+1
        tmpM['sample'] = np.repeat(range(int(plateW/2)),plateH*2)+1
        tmpM.loc[tmpM['sample'].isin([1,12]),'sample'] ='DMSO'
        tmpM = tmpM[['row','col','sample','value']]

        dmsoCtl = tmpM[(tmpM['sample'].isin(['DMSO'])) & (~tmpM['row'].isin(['A','B','O','P']))]

        
        lwr_b = stdDevbound(dmsoCtl['value'])[0]
        hir_b = stdDevbound(dmsoCtl['value'])[1]
        dmsoCtlwoOut = dmsoCtl[(dmsoCtl['value']>lwr_b) & (dmsoCtl['value']<hir_b)]
        outlier_def = dmsoCtl[~(dmsoCtl['value']>lwr_b) & (dmsoCtl['value']<hir_b)]
        
        print('For the file:')
        print(fileName)
        print('The outlier range is: ', [int(lwr_b),int(hir_b)])
        print('')
        print('Outliers within the DMSO are listed below:')
        print(outlier_def)
        print('')
        print('Outliers are removed by default option can be set to keep them or change their definition')
        if outliierCut == 1:
            tmpMnoOut = tmpM.drop(outlier_def.index, axis = 0)
        print('')
        print('------------------------------------------------------')

        tmpMnoOut['norm'] = tmpMnoOut['value']/np.mean(dmsoCtlwoOut['value'])
        tmpMnoOut = tmpMnoOut.rename(columns={'value':'value_'+dattype,'norm':'norm_'+dattype})

        return tmpMnoOut



    else:
        tmp = pd.read_csv(fileName)

    print('ERROR in reading the file the number of sample extracted is incorrect')
    return tmp
    
def combineData(mPath):

    dattype = 'ffluc'
    fflname = getFile(list(set(glob.glob(mPath+os.sep+'*'+dattype+'*.xlsx'))-set(glob.glob(mPath+os.sep+'~*'+dattype+'*.xlsx')))[0], dattype = dattype)

    dattype = 'nanolic'
    hibname = getFile(list(set(glob.glob(mPath+os.sep+'*'+dattype+'*.xlsx'))-set(glob.glob(mPath+os.sep+'~*'+dattype+'*.xlsx')))[0], dattype = dattype)

    combo = pd.merge(fflname,hibname,on=['row','col','sample'])

    return combo

def getThePlot(mPath, style='scatter'):
    print('Generating plots.....')

    fig, ax = plt.subplots(2,2,figsize=([16, 16]))


    combo = combineData(mPath)
    ##################################
    ## Get the plate map
    ##################################
    myVmin = min(np.array(axesParam(combo, False)).flatten())
    myVmax = max(np.array(axesParam(combo, False)).flatten())
    ax[0][0].imshow(getPlateFormat(combo, 'norm_ffluc', mPath), interpolation = 'nearest', aspect='auto', vmin=myVmin, vmax=myVmax) # cmap= 'bwr')
    ax[0][1].imshow(getPlateFormat(combo, 'norm_nanolic', mPath), interpolation = 'nearest', aspect='auto', vmin=myVmin, vmax=myVmax) # cmap= 'bwr')
    
    for i in ax[0]:
        i.set_xlabel('plate columns')
        i.set_ylabel('plate rows')

    ##################################
    ## Get the scatter
    ##################################
    if style == 'scatter':
        sns.scatterplot(data=combo, x="norm_ffluc", y="norm_nanolic", hue="sample", alpha=0.9, palette='Paired', ax=ax[1][0])

        dmsoCol = matplotlib.cm.get_cmap('Paired')(0)
        dmsoCol = matplotlib.colors.rgb2hex(dmsoCol)

        sns.scatterplot(data=combo[combo['sample']=='DMSO'], x="norm_ffluc", y="norm_nanolic", alpha=0.9, color=dmsoCol,ax=ax[1][1])

    elif style == 'kde':
        sns.kdeplot(data=combo, x="norm_ffluc", y="norm_nanolic", hue="sample", alpha=0.9, palette='Paired', ax=ax[1][0])

        dmsoCol = matplotlib.cm.get_cmap('Paired')(0)
        dmsoCol = matplotlib.colors.rgb2hex(dmsoCol)

        sns.kdeplot(data=combo[combo['sample']=='DMSO'], x="norm_ffluc", y="norm_nanolic", alpha=0.9, color=dmsoCol,ax=ax[1][1])

    xMark = stdDevbound(combo['norm_ffluc'])
    yMark = stdDevbound(combo['norm_nanolic'])

    for i in ax[1]:
        i.set_xlim(axesParam(combo)[0])
        i.set_ylim(axesParam(combo)[1])
        for j in [0,1]:
            # print(xMark[j])
            i.axvline(xMark[j], linestyle='dashed', color='red')
            i.axhline(yMark[j], linestyle='dashed', color='red')

    ax[0][0].set_title('Plate map of ffluc', fontsize= 12)
    ax[0][1].set_title('Plate map of nanolic', fontsize= 12)
    ax[1][0].set_title('Scatter for all the samples (red lines: 3*SD)', fontsize= 12)
    ax[1][1].set_title('Scatter for DMSO only (red lines: 3*SD)', fontsize= 12)

    saveName = mPath+os.sep+'figSummary_'+style
    plt.savefig(saveName+'.png')
    plt.savefig(saveName+'.pdf')

    ##################################
    ## Get the plot for individual samples
    ##################################
    
    fig, ax = plt.subplots(4,3,figsize=([16, 16]))

    for idx, (i, j) in enumerate(zip(ax.flatten(), combo['sample'].unique())):
        # print(idx, i, j)
        dmsoCol = matplotlib.cm.get_cmap('Paired')(idx)
        dmsoCol = matplotlib.colors.rgb2hex(dmsoCol)
        sns.scatterplot(data=combo[combo['sample']==j], x="norm_ffluc", y="norm_nanolic", alpha=0.9, color=dmsoCol,ax=i)
        sns.kdeplot(data=combo[combo['sample']==j], x="norm_ffluc", y="norm_nanolic", alpha=0.3, color=dmsoCol,ax=i)
        i.set_xlim(axesParam(combo)[0])
        i.set_ylim(axesParam(combo)[1])
        for k in [0,1]:
            # print(xMark[j])
            i.axvline(xMark[k], linestyle='dashed', color='red')
            i.axhline(yMark[k], linestyle='dashed', color='red')
            # i.set_title('Sample: '+str(j), fontsize= 12)
            i.text(0.1, 0.2, 'Sample: '+str(j), fontsize= 12)

    saveName = mPath+os.sep+'figSummaryInd_'+style
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
    sns.scatterplot(data=combo, x="norm_ffluc", y="norm_nanolic", s=80, color='black', alpha=0.2)
    sns.scatterplot(data=tmp, x="norm_ffluc", y="norm_nanolic", s=80, color='red', alpha=0.8)
    for i, j in tmp.iterrows():
        mytext = str(j['sample'])+': '+j['row']+str(j['col'])
        ax.text(j['norm_ffluc'], j['norm_nanolic'], mytext, fontsize = 12)

    print('')
    print('Table for data falling out of the 3*SD limits')
    print(tmp)

    tmp.to_csv(mPath+os.sep+'outerSample.csv')
    saveName = mPath+os.sep+'figSummaryOuterSample_'
    plt.savefig(saveName+'.png')
    plt.savefig(saveName+'.pdf')

def axesParam(data, atZero=True):
    '''
    data: correpsond to the combine data file 
    '''
    
    myXlim = [min(data['norm_ffluc'])*0.9, max(data['norm_ffluc'])*1.1]
    myYlim = [min(data['norm_nanolic'])*0.9, max(data['norm_nanolic'])*1.1]

    if atZero == True:
        myXlim = [0,1.1*max(max(data['norm_nanolic']),max(data['norm_ffluc']))]
        myYlim = [0,1.1*max(max(data['norm_nanolic']),max(data['norm_ffluc']))]

    return myXlim, myYlim

def stdDevbound(dataArray):
    '''
    dataArray: correspond to the array for which we want to obtain the limit of interest
    '''
    outlierCutoffValue = 3*np.std(dataArray)
    lwr_b = np.mean(dataArray)-outlierCutoffValue
    hir_b = np.mean(dataArray)+outlierCutoffValue

    return [lwr_b, hir_b]

def getPlateFormat(comboData, myVal, mPath):
    tmp = pd.pivot_table(comboData, values=myVal, index=['row'], columns=['col'])
    tmp.to_csv(mPath+os.sep+'normPlate_'+myVal+'.csv')
    tmp = tmp.to_numpy()

    return tmp

def getListOutofLim(combo):
    xMark = stdDevbound(combo['norm_ffluc'])
    yMark = stdDevbound(combo['norm_nanolic'])
    tmp = combo[~((combo['norm_ffluc']>xMark[0]) & (combo['norm_ffluc']<xMark[1]))]


    xmask = (combo['norm_ffluc']>xMark[0]) & (combo['norm_ffluc']<xMark[1])
    ymask = (combo['norm_nanolic']>yMark[0]) & (combo['norm_nanolic']<yMark[1])
    tmp = combo[~(xmask & ymask)]

    return tmp



    tmp = tmp[~((tmp['norm_nanolic']>yMark[0]) & (tmp['norm_ffluc']<yMark[1]))]
    combo[combo['norm_ffluc']]
outliierCut = 1

mPath = sys.argv[1]

# Check if path exits
if os.path.exists(mPath):
    print ("Folder exist")

# mPath = input('Drag the folder containing the file to analyze:')
# print(os.path.exists(mPath))
print(mPath)
# print(glob.glob(mPath+os.sep+'*.*'))
combo = combineData(mPath)
combo.to_csv(mPath+os.sep+'combinedOutput.csv')
getThePlot(mPath)

# 'C:/Users/Windows/Desktop/CamiloTest'