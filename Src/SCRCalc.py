import os
import numpy as np
import json as js
import pylab as pl
from scipy import integrate as Itg

# path='H:\Docs\Research\TCH\Second malignancy\BrstStudy'
# file1=path+'\\'+'1925520_Study2IMRT_DVH_1.txt'
# file2=path+'\\'+'1925520_Study2VMAT_DVH_1.txt'

def ReadDVHData(fp):
    dvhInfo={}
    fdata = open(fp, 'r')
    lines = fdata.readlines()
    fdata.close()
    nlines = np.size(lines)
    dvhInfo['PatID']=str(int(str(lines[0]).split('|')[0].split('~')[1]))
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

def WriteToJSON(data,fp):
    filepath=fp+'\\'+str(data['PatID'])+'.json'
    with open(filepath,'w') as OutFile:
        js.dump(data, OutFile)
        #print(filepath,':Written')

#Converts JSON to python dict
def JSONtoDict(JSONfile):
    JSONData=open(JSONfile)
    DVHData=js.load(JSONData)
    JSONData.close()
    return DVHData


def CalcOEDLinear(DoseBins,VolBins):
    totalVol=np.sum(VolBins)
    OED=[]
    # For linear DRC RED=D
    for x in range(0,np.size(VolBins),1):
        # Convert doses from cGy to GY
        RED=DoseBins[x]/100.0
        OED.append(VolBins[x]*RED)
    return np.sum(OED)*(1.0/totalVol)

def CalcOEDBell(DoseBins,VolBins,alpha):
    totalVol=np.sum(VolBins)
    OED=[]
    # For linear DRC RED=D
    for x in range(0,np.size(VolBins),1):
        # Convert doses from cGy to GY
        DoseBin=DoseBins[x]/100.0
        RED=DoseBin*(np.exp(-(alpha*DoseBin)))
        OED.append(VolBins[x]*RED)
    return np.sum(OED)*(1.0/totalVol)

def CalcOEDPlateau(DoseBins,VolBins,alpha):
    totalVol=np.sum(VolBins)
    OED=[]
    # For linear DRC RED=D
    for x in range(0,np.size(VolBins),1):
        # Convert doses from cGy to GY
        DoseBin=DoseBins[x]/100.0
        RED=(1.0-np.exp(-(alpha*DoseBin)))/(alpha)
        OED.append(VolBins[x]*RED)
    return np.sum(OED)*(1.0/totalVol)



def CalcOEDMechanistic(DoseBins,VolBins,alpha,ab,R,n):
    #RED=Risk equivalent dose
    totalVol = np.sum(VolBins)
    OED =[]
    for x in range(0, np.size(VolBins),1):
        if VolBins[x]>0.0:
            Dose=DoseBins[x]/100.0
            beta=alpha/ab
            alphaFr=alpha+(beta*(Dose/n))
            RED1=(np.exp(-alphaFr*Dose))/(alphaFr*R)
            RED2=(1.0-(2.0*R)+((R**2.0)*np.exp(alphaFr*Dose))-((1.0-R)**2.0)*np.exp(-(((alphaFr*R)*Dose)/(1.0-R))))
            #print(RED1,RED2,alphaFr)
            RED=RED1*RED2
            OED.append(VolBins[x]*RED)
    return np.sum(OED)*(1.0/totalVol)


def CalcMFAge(ageX,ageA,gammaE,gammaA):
    mu=np.exp(gammaE*(ageX-30)+(gammaA*np.log(ageA/70.0)))
    return mu


def CalcEAR(DoseBins,VolBins,alpha,ab,R,n,betaEAR,mu,genderCF):
    #RED=Risk equivalent dose
    totalVol = np.sum(VolBins)
    OED =[]
    for x in range(0, np.size(VolBins),1):
        if VolBins[x]>0.0:
            Dose=DoseBins[x]/100.0
            beta=alpha/ab
            alphaFr=alpha+(beta*(Dose/n))
            RED1=(np.exp(-alphaFr*Dose))/(alphaFr*R)
            RED2=(1-(2*R)+((R**2)*np.exp(alphaFr*Dose))-((1.0-R)**2)*np.exp(-(((alphaFr*R)*Dose)/(1-R))))
            #print(RED1,RED2,alphaFr)
            RED=RED1*RED2*betaEAR*mu
            OED.append(VolBins[x]*RED)
    return np.sum(OED)*(1.0/totalVol)*genderCF

def CalcEAR2(age,DoseBins,VolBins,alpha,ab,R,n,betaEAR,genderCF):
    #RED=Risk equivalent dose
    totalVol = np.sum(VolBins)
    mu=CalcMFAge(age,age,0.002,4.23)
    OED =[]
    Sax=CalcProbSurvival(age)
    for x in range(0, np.size(VolBins),1):
        if VolBins[x]>0.0:
            Dose=DoseBins[x]/100.0
            beta=alpha/ab
            alphaFr=alpha+(beta*(Dose/n))
            RED1=(np.exp(-alphaFr*Dose))/(alphaFr*R)
            RED2=(1-(2*R)+((R**2)*np.exp(alphaFr*Dose))-((1.0-R)**2)*np.exp(-(((alphaFr*R)*Dose)/(1-R))))
            #print(RED1,RED2,alphaFr)
            RED=RED1*RED2*betaEAR*mu
            OED.append(VolBins[x]*RED)
    EAR=((np.sum(OED)*(1.0/totalVol))*Sax*genderCF)
    if EAR<0.0:
        EAR=0
    return EAR


#Sax=S(age_a)/S(age_x)
def CalcLAR(DoseBins,VolBins,alpha,ab,R,n,betaEAR,genderCF,age1,age2):
    LAR=Itg.quad(CalcEAR2,age1,age2,args=(DoseBins,VolBins,alpha,ab,R,n,betaEAR,genderCF))
    return LAR[0]


def CalcgEUD(DoseBins,VolBins,n):
    totalVol = np.sum(VolBins)
    # print totalVol,'TV'
    if n == 0.0:
        n = 1e-7
    EUD = 0.0
    totalVol = np.sum(VolBins)
    # print totalVol,'TotalVol'
    for x in range(0, np.size(DoseBins), 1):
        if VolBins[x] > 0.0:
            EUD = EUD + ((DoseBins[x] ** (1.0 / n)) * (VolBins[x] / totalVol))
            # print EUD
    return (EUD**n)


#Calculates integral dose in litre-Gy units
def CalcIntegralDose(DoseBins,VolBins):
    MeanDose=CalcgEUD(DoseBins,VolBins,1)/100.0#(convert cGy to Gy)
    IntegralDose=MeanDose*(np.sum(VolBins)*0.001)
    return IntegralDose

#age exposed=39, polynomial fit only valid till age attained=80
def CalcProbSurvival(age_a):
    #prob=(-0.57173E-05*(age_a**3))+(0.6917706E-03*(age_a**2))-(0.0308322605*age_a)+1.4936202871
    prob=(6E-07*(age_a**3))-(5E-05*(age_a**2))-(0.0112*age_a)+1.0007
    #print(prob, "...")
    return prob



# #plot of probability of surviving till age x from age 39
# S_ae=[]
# age_e=[]
# for x in np.arange(39.0,81.0,1):
#     S_ae.append(CalcProbSurvival(x))
#     age_e.append(x)
#
# pl.plot(age_e,S_ae)
# pl.show()



# xLst=[]
# yLst1=[]
# yLst2=[]
# mu1=CalcMFAge(30,70,-0.056,6.9)
# mu2=CalcMFAge(30,70,-0.024,2.38)
# Vol=[2000.0]
# for x in np.arange(10,7000,1):
#     Dose=[x]
#     EAR1=CalcEAR(Dose,Vol,0.033,3.0,0.56,25,0.73,mu1)
#     EAR2 = CalcEAR(Dose, Vol, 0.219, 3.0, 0.06, 25, 3.8, mu2)
#     xLst.append(x/100.0)
#     yLst1.append(EAR1)
#     yLst2.append(EAR2)
#
# pl.plot(xLst,yLst1,linewidth=2.0)
# pl.hold(True)
# pl.plot(xLst,yLst2,linewidth=2.0)
# pl.legend(['Rectum','Bladder'],loc=2)
# pl.show()

# #Read and pack files
# folders=os.listdir(path)
# print("No. of folders found:",np.size(folders))
# for x in range(0,np.size(folders),1):
#     files=os.listdir(path+"\\"+folders[x])
#     curFile=path+"\\"+folders[x]+"\\"+files[1]#1 for vmat & 0 for IMRT
#     data=ReadDVHData(curFile)
#     print("Writing file...:",files[1],x,str(data['PatID']))
#     #WriteToJSON(data,"H:\Docs\Research\TCH\Second malignancy\JSONs\IMRT")
#     WriteToJSON(data,"H:\Desktop\BrstStudy-Updated Cord")


AgeList=[53,45,63,56,71,57,65,64,51,62,60,51,44,71,50,54,41,70,57,53,60,54,41,51,83,69,59,43,72,54,62,40,52,66,43,48,51,65,48,58,54,42,51,49,43,53,72,42,54,56]
#
# path="H:\Docs\Research\TCH\Second malignancy\JSONs\VMAT"
# files=os.listdir(path)
# for x in range(0,np.size(files),1):
#     curFile=path+"\\"+files[x]
#     DVHs=JSONtoDict(curFile)
#     try:
#         DVH1=DVHs['Contralateral Breast']
#         doses=DVH1['DoseBins']
#         Vols=DVH1['VolBins']
#         #OED=CalcOEDLinear(doses,Vols)
#         #OED=CalcOEDPlateau(doses,Vols,0.009)
#         #OED=CalcOEDBell(doses,Vols,0.009)
#         #OED=CalcOEDMechanistic(doses,Vols,0.018,3.0,0.93,25)
#         #OED=CalcIntegralDose(doses,Vols)
#         mu=CalcMFAge(30.0,70.0,-0.037,1.7)
#         EAR=CalcEAR(doses,Vols,0.044 ,3.0,0.15,25,8.0,mu)
#         print(EAR)
#     except:
#         print("ROI not found:",curFile)


# curFile="H:\\Docs\\Research\\TCH\\Second malignancy\\JSONs\\VMAT\\2003750.json"
# DVHs=JSONtoDict(curFile)
# DVH1=DVHs['Contralateral Breast']
# doses=DVH1['DoseBins']
# vols=DVH1['VolBins']
# pl.plot(doses,vols)
# pl.show()



# #Local dose-effect curves for various models
# doses=np.arange(0,10000,50)
# vols=doses*0.0+5
# OEDLinear=[]
# OEDBell=[]
# OEDPlat=[]
# OEDMechHalfRepop=[]
# OEDMechNoRepop=[]
# OEDMechFullRepop=[]
#
# for x in range(0,np.size(doses),1):
#     OEDLinear.append(CalcOEDLinear([doses[x]],[5]))
#     OEDBell.append(CalcOEDBell([doses[x]],[5],0.041))
#     OEDPlat.append(CalcOEDPlateau([doses[x]],[5],0.041))
#     OEDMechHalfRepop.append(CalcOEDMechanistic([doses[x]],[5],0.044,10,0.5,25))
#     OEDMechNoRepop.append(CalcOEDMechanistic([doses[x]], [5], 0.044, 10, 0.01,25))
#     OEDMechFullRepop.append(CalcOEDMechanistic([doses[x]], [5], 0.044, 10, 1.0,25))


# pl.plot(doses/100.0,OEDLinear,linewidth=4)
# pl.hold(True)
# pl.plot(doses/100.0,OEDBell,linewidth=4,linestyle='-',marker='o')
# pl.hold(True)
# pl.plot(doses/100.0,OEDPlat,linewidth=4,linestyle='-',marker='+')
# pl.hold(True)
# pl.plot(doses/100,OEDMechHalfRepop,linewidth=4,linestyle='-.')
# pl.hold(True)
# pl.plot(doses/100,OEDMechNoRepop,linewidth=4,linestyle='--')
# pl.hold(True)
# pl.plot(doses/100,OEDMechFullRepop,linewidth=4,linestyle=':')
# pl.legend(['Linear','Bell','Plateau','Mechanistic (R=0.5)','Mechanistic (R=0.01)','Mechanistic (R=1.0)'])
# pl.grid(False)
# pl.xticks(fontsize=15)
# pl.yticks(fontsize=15)
# pl.xlabel('Dose (Gy)',fontsize=20)
# pl.ylabel('OED (Gy)',fontsize=20)
# pl.show()
# #End of loacl dose-effect curves for various models


#Calculate LAR #
path="H:\Docs\Research\TCH\Second malignancy\JSONs\IMRT"
files=os.listdir(path)
LARList=[]
RRList=[]
for x in range(0,np.size(files),1):
    curFile=path+"\\"+files[x]
    DVHs=JSONtoDict(curFile)
    DVH1=DVHs['Rt Lung']
    doses=DVH1['DoseBins']
    Vols=DVH1['VolBins']
    FinalLAR=0
    LAR=CalcLAR(doses,Vols,0.042,3.0,0.83,25,8.0,1.17,AgeList[x],75)
    if LAR<0.0:
        LAR=0.0
    print(LAR)
    RR=(LAR+49.4)/49.4
    RRList.append(RR)
    LARList.append(LAR)
# print(np.mean(LARList),np.std(LARList))
# print(np.mean(RRList),np.std(RRList))










