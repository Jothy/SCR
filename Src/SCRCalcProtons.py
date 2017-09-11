import os
import numpy as np
import numpy.random as npr
import pylab
import json as js
import pylab as pl
from scipy import integrate as Itg


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


def BootstrapCI(data, num_samples, statistic, alpha):
    """Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic."""
    samples = npr.choice(data, size=(num_samples, len(data)), replace=True)
    stat = np.sort(statistic(samples, 1))
    return (stat[int((alpha/2.0)*num_samples)],
            stat[int((1-alpha/2.0)*num_samples)])


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

#Converts JSON to python dict
def JSONtoDict(JSONfile):
    JSONData=open(JSONfile)
    DVHData=js.load(JSONData)
    JSONData.close()
    return DVHData



# data1=[16.69,	19.91,	19.82,	18.00,	17.07,	14.82,	17.88,	19.30,	15.63,	18.04]
# data2=[13.00,	14.19,	13.76,	11.75,	11.60,	9.10,	13.85,	13.92,	9.88,	9.90]
# data3=[12.75,	14.22,	13.77,	11.78,	11.22,	9.39,	12.38,	13.67,	10.57,	9.07]
# data4=[9.89,	13.30,	11.07,	13.32,	9.60,	9.63,	9.75,	11.22,	12.15,	9.74]
# data5=[0.64, 1.06,	1.09,	0.36,	0.33,	0.90,	0.23,	1.24,	1.06,	0.32]
#
#
# low,high=BootstrapCI(data1,1000000,np.mean,0.05)
# print(np.round(low,2),np.round(high,2), "95% CI")
#
# low,high=BootstrapCI(data2,1000000,np.mean,0.05)
# print(np.round(low,2),np.round(high,2), "95% CI")
#
# low,high=BootstrapCI(data3,1000000,np.mean,0.05)
# print(np.round(low,2),np.round(high,2), "95% CI")
#
# low,high=BootstrapCI(data4,1000000,np.mean,0.05)
# print(np.round(low,2),np.round(high,2), "95% CI")
#
# low,high=BootstrapCI(data5,1000000,np.mean,0.05)
# print(np.round(low,2),np.round(high,2), "95% CI")

path="D:\\Projects\\SCR\\Data\\CSI_JSONs\\Protons"
files=os.listdir(path)
for x in range(0,np.size(files),1):
     curFile=path+"\\"+files[x]
     #print(curFile)
     DVHs=JSONtoDict(curFile)
     try:
         DVH1=DVHs['Bowel Bag ']
         doses=DVH1['DoseBins']
         Vols=DVH1['VolBins']
         #OED=CalcOEDLinear(doses,Vols)
         #MeanDose=CalcgEUD(doses,Vols,1)
#         #OED=CalcOEDPlateau(doses,Vols,0.009)
#         #OED=CalcOEDBell(doses,Vols,0.009)
         OED=CalcOEDMechanistic(doses,Vols,0.089,3.0,0.17,20)
#         #OED=CalcIntegralDose(doses,Vols)
#         mu=CalcMFAge(30.0,70.0,-0.037,1.7)
#         EAR=CalcEAR(doses,Vols,0.044 ,3.0,0.15,25,8.0,mu)
         print(OED)
     except:
         print("ROI not found:",curFile)
print('************************************')
#
# path="D:\\Projects\\SCR\\Data\\CSI_JSONs\\IMRT"
# files=os.listdir(path)
# for x in range(0,np.size(files),1):
#      curFile=path+"\\"+files[x]
#      #print(curFile)
#      DVHs=JSONtoDict(curFile)
#      try:
#          DVH1=DVHs['Stomach']
#          doses=DVH1['DoseBins']
#          Vols=DVH1['VolBins']
#          #OED=CalcOEDLinear(doses,Vols)
#          MeanDose=CalcgEUD(doses,Vols,1)
# #         #OED=CalcOEDPlateau(doses,Vols,0.009)
# #         #OED=CalcOEDBell(doses,Vols,0.009)
# #         #OED=CalcOEDMechanistic(doses,Vols,0.018,3.0,0.93,25)
# #         #OED=CalcIntegralDose(doses,Vols)
# #         mu=CalcMFAge(30.0,70.0,-0.037,1.7)
# #         EAR=CalcEAR(doses,Vols,0.044 ,3.0,0.15,25,8.0,mu)
#          print(MeanDose/100.0)
#      except:
#          print("ROI not found:",curFile)
# print('************************************')
#
# path="D:\\Projects\\SCR\\Data\\CSI_JSONs\\VMAT"
# files=os.listdir(path)
# for x in range(0,np.size(files),1):
#      curFile=path+"\\"+files[x]
#      #print(curFile)
#      DVHs=JSONtoDict(curFile)
#      try:
#          DVH1=DVHs['Stomach']
#          doses=DVH1['DoseBins']
#          Vols=DVH1['VolBins']
#          #OED=CalcOEDLinear(doses,Vols)
#          MeanDose=CalcgEUD(doses,Vols,1)
# #         #OED=CalcOEDPlateau(doses,Vols,0.009)
# #         #OED=CalcOEDBell(doses,Vols,0.009)
# #         #OED=CalcOEDMechanistic(doses,Vols,0.018,3.0,0.93,25)
# #         #OED=CalcIntegralDose(doses,Vols)
# #         mu=CalcMFAge(30.0,70.0,-0.037,1.7)
# #         EAR=CalcEAR(doses,Vols,0.044 ,3.0,0.15,25,8.0,mu)
#          print(MeanDose/100.0)
#      except:
#          print("ROI not found:",curFile)
# print('************************************')
#
# path="D:\\Projects\\SCR\\Data\\CSI_JSONs\\Tomo"
# files=os.listdir(path)
# for x in range(0,np.size(files),1):
#      curFile=path+"\\"+files[x]
#      #print(curFile)
#      DVHs=JSONtoDict(curFile)
#      try:
#          DVH1=DVHs['Stomach ']
#          doses=DVH1['DoseBins']
#          Vols=DVH1['VolBins']
#          #OED=CalcOEDLinear(doses,Vols)
#          MeanDose=CalcgEUD(doses,Vols,1)
# #         #OED=CalcOEDPlateau(doses,Vols,0.009)
# #         #OED=CalcOEDBell(doses,Vols,0.009)
# #         #OED=CalcOEDMechanistic(doses,Vols,0.018,3.0,0.93,25)
# #         #OED=CalcIntegralDose(doses,Vols)
# #         mu=CalcMFAge(30.0,70.0,-0.037,1.7)
# #         EAR=CalcEAR(doses,Vols,0.044 ,3.0,0.15,25,8.0,mu)
#          print(MeanDose/100.0)
#      except:
#          print("ROI not found:",curFile)
# print('************************************')
#
#
# path="D:\\Projects\\SCR\\Data\\CSI_JSONs\\Protons"
# files=os.listdir(path)
# for x in range(0,np.size(files),1):
#      curFile=path+"\\"+files[x]
#      #print(curFile)
#      DVHs=JSONtoDict(curFile)
#      try:
#          DVH1=DVHs['Stomach ']
#          doses=DVH1['DoseBins']
#          Vols=DVH1['VolBins']
#          #OED=CalcOEDLinear(doses,Vols)
#          MeanDose=CalcgEUD(doses,Vols,1)
# #         #OED=CalcOEDPlateau(doses,Vols,0.009)
# #         #OED=CalcOEDBell(doses,Vols,0.009)
# #         #OED=CalcOEDMechanistic(doses,Vols,0.018,3.0,0.93,25)
# #         #OED=CalcIntegralDose(doses,Vols)
# #         mu=CalcMFAge(30.0,70.0,-0.037,1.7)
# #         EAR=CalcEAR(doses,Vols,0.044 ,3.0,0.15,25,8.0,mu)
#          print(MeanDose/100.0)
#      except:
#          print("ROI not found:",curFile)


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

