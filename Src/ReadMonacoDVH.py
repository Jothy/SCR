import os
import numpy as np
import json as js
import pylab as pl
from scipy import integrate as Itg
import xlrd


def ReadDVHDataMonaco(fp):
    dvhInfo={}
    fdata = open(fp, 'r')
    lines = fdata.readlines()
    fdata.close()
    nlines = np.size(lines)
    dvhInfo['PatID']=str(lines[0]).split('|')[0].split('~')[1]
    dvhInfo['PlanName']=str(lines[0]).split('|')[1].split(':')[1]
    dvhInfo['DoseUnits'] = str(lines[0]).split('|')[4].split(':')[1]
    dvhInfo['VolumeUnits'] = str(lines[0]).split('|')[5].split(':')[1][0:4]
    #print(dvhInfo['PatID'],dvhInfo['PlanName'],dvhInfo['DoseUnits'],dvhInfo['VolumeUnits'])
    #Remove headers & footers
    linesFiltered=lines[3:(nlines-3)]
    doses=[]
    volumes=[]
    for x in range(0,np.size(linesFiltered),1):
       doses.append(np.float(linesFiltered[x].split('                   ')[1]))
       volumes.append(np.float(linesFiltered[x].split('                   ')[2]))
    dvhStartStopIndices=[item for item in range(len(doses)) if doses[item] == 0]
    dvhStartStopIndices.append(np.size(linesFiltered))
    for x in range(0,np.size(dvhStartStopIndices)-1,1):
        curDVH = {}
        curDVH['DoseBins']=doses[dvhStartStopIndices[x]:dvhStartStopIndices[x+1]]
        curDVH['VolBins']=volumes[dvhStartStopIndices[x]:dvhStartStopIndices[x+1]]
        ROI=linesFiltered[dvhStartStopIndices[x]].split('                   ')[0]
        dvhInfo[str(ROI)]=curDVH
        #print(np.sum(curDVH['VolBins']),ROI)
    return dvhInfo

def ReadTomo(fp):
    DVHData = list()
    curFile = xlrd.open_workbook(fp)
    curSheet = curFile.sheet_by_name(pageName)
    FrsCol = curSheet.col_values(FrsColNum)  # Frs column
    MLDCol = curSheet.col_values(MLDColNum)  # MLD column
    RPCol = curSheet.col_values(RPColNum)  # RP column
    PDCol = curSheet.col_values(FrsColNum - 1)  # Prescription dose


def WriteToJSON(data,fp):
    filepath=fp+'\\'+str(data['PatID'])+'_VMAT'+'.json'
    with open(filepath,'w') as OutFile:
        js.dump(data, OutFile)
        #print(filepath,':Written')

#Converts JSON to python dict
def JSONtoDict(JSONfile):
    JSONData=open(JSONfile)
    DVHData=js.load(JSONData)
    JSONData.close()
    return DVHData


# fp="D:\Projects\SCR\Data\CSI10\CSI10_VMAT_DVH_1.txt"
# fw="D:\Projects\SCR\Data\CSI10"
# data=ReadDVHDataMonaco(fp)
# WriteToJSON(data,fw)

curFile=xlrd.open_workbook('D:\Projects\SCR\Data\CSI1\CSI1_Tomo.xlsx')
curSheet=curFile.sheet_by_name('CSI1_Tomo')










