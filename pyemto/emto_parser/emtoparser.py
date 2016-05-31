import matplotlib.pyplot as plt
from os import popen
import numpy as np
import pandas as pd

#pd.options.display.max_columns = 5200
#pd.set_option('notebook_repr_html',True)
#pd.set_option('display.width', 1000)
#pd.set_option('display.max_colwidth', 100)

All = slice(None)

#Kb = 8.617332478e-5*0.073498618 #Ry/K
Kb = 1.3806488e-23*1000#mJ/K
Bohr2Anstr = 0.529177249
Anstr2M = 1e-10
Ry2mJ = 2.1798741e-18 * 1000
RyPbohr32Gpa = 14710.5
RyPbohr22mJPm2 = Ry2mJ/(Bohr2Anstr*Anstr2M)**2

nameparserexample = r'KFCD-(?P<STR>...).*Ni(?P<NiCC>[0-8]+)-(?P<SWS0>[0-9]\.[0-9]+)-.*-(?P<REX>[0-9]+)-.*'

def Sws2T(sws,sws0,alpha):
    return (sws-sws0)/sws0/alpha            

def str2num(string):
    '''
    change string to int, float, complex or itself
    
    example:
        >>> stringlist = ["1",'Fe3','4.87','1e-7','4.6+9j']
        >>> map(ok, stringlist)
        [1, 'Fe3', 4.87, 1e-07, (4.6+9j)]
    '''
    if not isinstance(string, basestring):
        #raise ValueError('input must be string, yours is a '+type(string).__name__)
        return string
    try:
        return int(string)
    except ValueError:
        try:
            return float(string)
        except ValueError:
            try:
                return complex(string)
            except ValueError:            
                return string

class EMTOPARSER:

    def __init__(self,filename,xc="PBE",suffix="kfcd"):
        """
        parse the kfcd output
        after init you can get a pd.DataFrame self.EN
        usage:
            
            FeCrNi = EMTOPARSER("./*","QNA","dat")
                you can access DataFrame: FeCrNi.EN
                FeCrNi.EN have FN Etot Esite and SWS as columns
           
            if they are Magnetic system you can call:
            
            FeCrNi.MagTables()
                you can access DataFrame: FeCrNi.ES and FeCrNi.CM
                FeCrNi.EN have FN Etot Esite SWS and S as columns
                FeCrNi.CM have Mag and CC for each IQ ITA and element
        """
        self.filename = filename + '.' + suffix
        self.xc = xc
        self.EnergyColume = [0, 2, 4, 7]
        self.EnergyColumeName = ["FN","Etot","Esite","SWS"]
        self.StructureColume = [0,1]
        self.StructureColumeName = ["FN","str"]
        self.MagneticColume = [0, 6, 8]
        self.MagneticColumeName = ["FN","IQ","Mag"]
        self.ConcentrationColume = [0, 3, 6, 9, 10]
        self.ConcentrationColumeName = ["FN","IQ","ITA","CC",'ELE']
        self.nameparser = None
        self.EN = self.Df(self.Energy(),self.EnergyColume,self.EnergyColumeName)
        self.STR = self.Df(self.Structure(),self.StructureColume,self.StructureColumeName)
        self.EN = self.EN.join(self.STR.set_index(["FN"]),on=["FN"]).applymap(str2num)

    def Energy(self):
        cmd =  "grep -H TOT-{} {}".format(self.xc,self.filename) 
        return [i.split() for i in popen(cmd).readlines()]

    def Structure(self):
        cmd = r'grep -H mdl {} |sed "s: .*/\([^/].*\).mdl: \1:" '.format(self.filename)
        return [i.split() for i in popen(cmd).readlines()]

    def Mag(self):
        cmd =  "grep -H Mag {}".format(self.filename) 
        return [i.split() for i in popen(cmd).readlines()]

    def Concentration(self):
        cmd =  'grep -H \'IQ.*)\' {}'.format(self.filename) 
        return [i.split() for i in popen(cmd).readlines()]

    def Df(self,list,col,colname):
        return pd.DataFrame([[i[x] for x in col] for i in list],columns=colname)
         
    def MagTables(self):
        mag_df = self.Df(self.Mag(),self.MagneticColume,self.MagneticColumeName)
        #self.STR.reset_index(inplace=True)
        #self.STR.set_index(["FN"],inplace=True)
        #mag_df = mag_df.join(self.STR,on=["FN"])
        #
        self.EN.reset_index(inplace=True)
        self.EN.set_index(["FN"],inplace=True)
        mag_df = mag_df.join(self.EN,on=["FN"])
        if self.nameparser is None:
            pass
        else:
            fn_df = pd.DataFrame(mag_df.FN.str.extract(self.nameparser))
            fnc = fn_df.columns.values.tolist()
            fnc.append("FN")
            fnc.append("str")
            mag_df = mag_df.join(fn_df)
        
        cn_df  = self.Df(self.Concentration(),self.ConcentrationColume,self.ConcentrationColumeName)
        cn_df.ELE = cn_df.ELE.str.replace("(","")
        #self.cn = cn_df
        #self.mag_df = mag_df
        cm_df  = cn_df.join(mag_df,rsuffix='r')
        if (cm_df.FN != cm_df.FNr).any():
            print "big problem !!!"
            raise ValueError
        cm_df.drop(["FNr","IQr"],axis=1,inplace=True)

        cm_df = pd.pivot_table(cm_df.applymap(str2num),index=["FN","str","SWS","Etot"],values=["Mag","CC"],columns=["IQ","ITA","ELE"])
        #may be a bug here, we can not get ELE any more, but it seem no harm.

        #if self.nameparser is None:
        #    cm_df = pd.pivot_table(cm_df.applymap(str2num),index=["FN","str","SWS","Etot"],values=["Mag","CC","ELE"],columns=["IQ","ITA"])
        #else:
        #    cm_df = pd.pivot_table(cm_df.applymap(str2num),index=fnc,values=["Mag","CC","ELE"],columns=["IQ","ITA"])
            
        self.CM = cm_df.reset_index()
        
        #calculate Magnetic Entropy
        #S = np.sum(np.log(np.abs(cm_df.Mag)+1.)*cm_df.CC,axis=1)
        #self.S_df = pd.DataFrame(S,columns=["S"])
        #self.S_df.reset_index(inplace=True)
        #self.S_df.set_index(["FN","str"],inplace=True)
        #self.ES = self.EN.join(self.S_df,on=["FN","str"])


if __name__=='__main__' :
    test = EMTOPARSER("kfcd_output/al0.20_co0.20_cr0.20_fe0.20_ni0.20_bcc_2.550000",suffix='prn')
    print test.Energy()
    print test.Structure()
    print test.Concentration()
    test.MagTables()
    print test.CM
    #print test.ES
    #print test.EN
    #print test
    #test.Df().plot()
    #test.EN.plot()
    #plt.show()
