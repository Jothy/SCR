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



def ReadTomoDVHFromXLS(fp):
    dvhInfo = {}
    dvhInfo['PatID'] ='12345'
    dvhInfo['PlanName'] = 'CSI'
    dvhInfo['DoseUnits'] ='Gy'
    dvhInfo['VolumeUnits'] ='cc'
    curFile = xlrd.open_workbook(fp)
    curSheet=curFile.sheet_by_index(0)
    # Read the dose column first
    DoseCol = curSheet.col_values(0)
    DoseCol = DoseCol[2:np.size(DoseCol)]
    [float(i) for i in DoseCol]

    for x in range(1, curSheet.ncols, 1):
        curDVH = {}
        ROINum = x
        ROI = curSheet.col_values(ROINum)
        ROIName = ROI[1].split(sep='(')[0]
        #print(ROIName)
        TotalVol = np.float(ROI[1].split(sep='(')[2].split(sep=':')[1].split(sep=')')[0])
        # print(TotalVol)
        VolCol = ROI[2:np.size(ROI)]
        while '' in VolCol:
            VolCol.remove('')
        [float(i) for i in VolCol]
        # Convert cumulative DVH to differential DVH
        DiffVols = []
        DiffDoses = []
        for y in range(0, np.size(VolCol) - 1, 1):
            DiffVols.append(((VolCol[y] - VolCol[y + 1]) * TotalVol)/100.0)
            DiffDoses.append((DoseCol[y] + DoseCol[y + 1])/2.0)

        curDVH['DoseBins'] =DiffDoses
        curDVH['VolBins'] = DiffVols
        dvhInfo[str(ROIName)] = curDVH
    return dvhInfo


def WriteToJSON(data,fp):
    filepath=fp+'\\'+str(data['PatID'])+'_Tomo'+'.json'
    with open(filepath,'w') as OutFile:
        js.dump(data, OutFile)
        #print(filepath,':Written')

#Converts JSON to python dict
def JSONtoDict(JSONfile):
    JSONData=open(JSONfile)
    DVHData=js.load(JSONData)
    JSONData.close()
    return DVHData



fp='D:\\Projects\\SCR\\Data\\CSI5_Tomo.xlsx'
Data=ReadTomoDVHFromXLS(fp)
WriteToJSON(Data,'D:\\Projects\\SCR\\Data')







# fp="D:\Projects\SCR\Data\CSI10\CSI10_VMAT_DVH_1.txt"
# fw="D:\Projects\SCR\Data\CSI10"
# data=ReadDVHDataMonaco(fp)
# WriteToJSON(data,fw)

# #plot diff DVH for testing
# pl.plot(DiffDoses1,DiffVols,linewidth=3.0)
# pl.xlim(0)
# pl.ylim(0)
# pl.show()


# print(ROIName,TotalVol)
# pl.plot(DoseCol[0:np.size(VolCol)],VolCol,linewidth=3.0)
# pl.xlim(0)
# pl.ylim(0,100)
# pl.grid(True)
# pl.show()














