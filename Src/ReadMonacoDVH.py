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

fp='D:\Projects\SCR\Data\CSI5\CSI5_Tomo.xlsx'
def ReadDVHDataTomoProton(fp):
    dvhInfo ={}

    curFile = xlrd.open_workbook(fp)
    sheetName1=fp.split('\\')[-1]
    sheetName2=sheetName1.split('.')[0]
    curSheet = curFile.sheet_by_name(sheetName2)

    PatID=sheetName2.split('_')[0]

    dvhInfo['PatID'] = PatID
    dvhInfo['PlanName'] = sheetName2
    dvhInfo['DoseUnits'] = "cGy"
    dvhInfo['VolumeUnits'] = "cc"

    # Read the dose column first
    DoseCol = curSheet.col_values(0)
    DoseCol = DoseCol[2:np.size(DoseCol)]
    [float(i) for i in DoseCol]  # convert to float
    DoseCol = [x * 100.0 for x in DoseCol]  # Convert from Gy to cGy
    #print(curSheet.ncols, ":Cols")
    for x in range(1, curSheet.ncols, 1):
        ROINum = x
        ROI = curSheet.col_values(ROINum)
        ROIName = ROI[1].split(sep='(')[0]
        TotalVol = np.float(ROI[1].split(sep='(')[2].split(sep=':')[1].split(sep=')')[0])
        VolCol = ROI[2:np.size(ROI)]
        while '' in VolCol:
            VolCol.remove('')
        [float(i) for i in VolCol]
        DiffDose = []
        DiffVols = []
        DiffDose, DiffVols = CumToDiffDVH(DoseCol, VolCol, TotalVol)
        pl.plot(DiffDose, DiffVols, linewidth=3.0)

        curDVH = {}
        curDVH['DoseBins'] = DiffDose
        curDVH['VolBins'] = DiffVols
        dvhInfo[str(ROIName)] =curDVH
    return dvhInfo


def WriteToJSON(data,fp):
    filepath=fp+'\\'+str(data['PatID'])+'_Proton'+'.json'
    with open(filepath,'w') as OutFile:
        js.dump(data, OutFile)
        #print(filepath,':Written')

#Converts JSON to python dict
def JSONtoDict(JSONfile):
    JSONData=open(JSONfile)
    DVHData=js.load(JSONData)
    JSONData.close()
    return DVHData

#Dose in cGY,VOls in %
def CumToDiffDVH(Doses,Vols,TotalVol):
    nbins=np.size(Vols)
    DiffVols=[]
    DiffDoses=[]
    for x in range(0,nbins-1,1):
        DiffVols.append(((Vols[x]-Vols[x+1])*TotalVol)/100.0)
        DiffDoses.append((Doses[x]+Doses[x+1])/2.0)
    return DiffDoses,DiffVols




# fp="D:\Projects\SCR\Data\CSI10\CSI10_VMAT_DVH_1.txt"
# fw="D:\Projects\SCR\Data\CSI10"
# data=ReadDVHDataMonaco(fp)
# WriteToJSON(data,fw)


fp='D:\Projects\SCR\Data\CSI1\CSI1_Proton.xlsx'
fw="D:\Projects\SCR\Data\CSI1"
data=ReadDVHDataTomoProton(fp)
WriteToJSON(data,fw)

