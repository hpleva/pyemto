import matplotlib.pyplot as plt
from os import popen
import numpy as np
import pandas as pd
import pymatgen as pmg

#pd.options.display.max_columns = 5200
#pd.set_option('notebook_repr_html',True)
#pd.set_option('display.width', 1000)
#pd.set_option('display.max_colwidth', 100)

All = slice(None)

#Kb = 8.617332478e-5*0.073498618 #Ry/K
#Kb = 1.3806488e-23*1000#mJ/K
kb = 8.6173324E-5 #eV/K
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

    def __init__(self,KGRN_filenames,KFCD_filenames,xc="PBE",suffix="prn",mag=False,DLM=False):
        """
        parse KGRN and KFCD output
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
        self.KGRN_filenames = KGRN_filenames + '.' + suffix        
        self.KFCD_filenames = KFCD_filenames + '.' + suffix
        self.xc = xc
        self.mag = mag
        self.DLM = DLM
        self.EnergyColumn = [0, 2, 4, 7]
        self.EnergyColumnName = ["FN","Etot","Eatom","SWS"]
        self.StructureColumn = [0,1]
        self.StructureColumnName = ["FN","str"]
        self.MagneticColumn = [0, 6, 8]
        self.MagneticColumnName = ["FN","IQ","Mag"]
        self.ConcentrationColumn = [0, 3, 6, 9, 10]
        self.ConcentrationColumnName = ["FN","IQ","ITA","Conc",'Elem']
        self.DOSColumn = [0]
        self.DOSColumnName = ['DOS']
        self.nameparser = None
        self.EN = self.Df(self.Energy(),self.EnergyColumn,self.EnergyColumnName)
        self.STR = self.Df(self.Structure(),self.StructureColumn,self.StructureColumnName)
        self.EN = self.EN.join(self.STR.set_index(["FN"]),on=["FN"]).applymap(str2num)

    def Energy(self):
        cmd =  "grep -H TOT-{} {}".format(self.xc,self.KFCD_filenames) 
        return [i.split() for i in popen(cmd).readlines()]

    def Structure(self):
        cmd = r'grep -H mdl {} |sed "s: .*/\([^/].*\).mdl: \1:" '.format(self.KFCD_filenames)
        return [i.split() for i in popen(cmd).readlines()]

    def Mag(self):
        cmd =  "grep -H Mag {}".format(self.KFCD_filenames) 
        return [i.split() for i in popen(cmd).readlines()]

    def Concentration(self):
        cmd =  'grep -H \'IQ.*)\' {}'.format(self.KFCD_filenames) 
        return [i.split() for i in popen(cmd).readlines()]

    def DOSEF(self):
        cmd =  'grep -H \'Dos(Ef)\' {}'.format(self.KGRN_filenames)
        cmd_output = popen(cmd).readlines()
        final_output = []
        if self.mag == False:
            for i in cmd_output:
                zero_line = True
                dos_line = i.split()[3:]
                for j in dos_line:
                    if j != '0.000000':
                        zero_line = False
                        break
                if zero_line == False:
                    final_output.append([i.split()[-1]])
        elif self.mag == True:
            for i in cmd_output:
                dos_line = i.split()
                if len(dos_line) == 4:
                    if dos_line[-1] != '0.000000':
                        final_output.append([dos_line[-1]])
        return final_output
        
    def Df(self,list,col,colname):
        return pd.DataFrame([[i[x] for x in col] for i in list],columns=colname)
         
    def create_df(self):
        self.EN.reset_index(inplace=True)
        self.EN.set_index(["FN"],inplace=True)
        #

        # Extract DOS(Ef) out of KGRN output files
        self.dos_df = self.Df(self.DOSEF(),self.DOSColumn,self.DOSColumnName)

        if self.mag == False:
                
            cm_df  = self.Df(self.Concentration(),self.ConcentrationColumn,self.ConcentrationColumnName)
            cm_df.Elem = cm_df.Elem.str.replace("(","")

            cm_df = cm_df.join(self.EN,on=["FN"])
            
            #cm_df  = cn_df.join(mag_df,rsuffix='r')
            #if (cm_df.FN != cm_df.FNr).any():
            #    print "big problem !!!"
            #    raise ValueError
            #cm_df.drop(["FNr","IQr"],axis=1,inplace=True)
            
            cm_df = pd.pivot_table(cm_df.applymap(str2num),index=["FN","str","SWS","Etot","Eatom"],values=["Conc"],columns=["IQ","ITA","Elem"])
            #may be a bug here, we can not get Elem any more, but it seem no harm.
            
            #if self.nameparser is None:
            #    cm_df = pd.pivot_table(cm_df.applymap(str2num),index=["FN","str","SWS","Etot"],values=["Mag","CC","Elem"],columns=["IQ","ITA"])
            #else:
            #    cm_df = pd.pivot_table(cm_df.applymap(str2num),index=fnc,values=["Mag","CC","Elem"],columns=["IQ","ITA"])
            
        elif self.mag == True:
            mag_df = self.Df(self.Mag(),self.MagneticColumn,self.MagneticColumnName)
            mag_df = mag_df.join(self.EN,on=["FN"])
        
            if self.nameparser is None:
                pass
            else:
                fn_df = pd.DataFrame(mag_df.FN.str.extract(self.nameparser))
                fnc = fn_df.columns.values.tolist()
                fnc.append("FN")
                fnc.append("str")
                mag_df = mag_df.join(fn_df)
        
            cn_df  = self.Df(self.Concentration(),self.ConcentrationColumn,self.ConcentrationColumnName)
            cn_df.Elem = cn_df.Elem.str.replace("(","")
            #self.cn = cn_df
            #self.mag_df = mag_df
            cm_df  = cn_df.join(mag_df,rsuffix='r')
            if (cm_df.FN != cm_df.FNr).any():
                print "big problem !!!"
                raise ValueError
            cm_df.drop(["FNr","IQr"],axis=1,inplace=True)
            
            cm_df = pd.pivot_table(cm_df.applymap(str2num),index=["FN","str","SWS","Etot","Eatom"],values=["Mag","Conc"],columns=["IQ","ITA","Elem"])
            #may be a bug here, we can not get Elem any more, but it seem no harm.

            #if self.nameparser is None:
            #    cm_df = pd.pivot_table(cm_df.applymap(str2num),index=["FN","str","SWS","Etot"],values=["Mag","CC","Elem"],columns=["IQ","ITA"])
            #else:
            #    cm_df = pd.pivot_table(cm_df.applymap(str2num),index=fnc,values=["Mag","CC","Elem"],columns=["IQ","ITA"])

        # Some trick
        self.CM = cm_df.reset_index()
        
        # Calculate configurational entropy
        if self.DLM == False:
            Sconf = -kb*np.sum(np.log(self.CM.Conc[1])*self.CM.Conc[1],axis=1)
            self.Sconf_df = pd.DataFrame(Sconf,columns=["Sconf"])
            #self.Sconf_df.reset_index(inplace=True)
            #self.Sconf_df = self.Sconf_df[["Sconf"]]
        elif self.DLM == True:
            Sconf = -kb*np.sum(np.log(2*self.CM.Conc[1].ix[:,1::2])*2*self.CM.Conc[1].ix[:,1::2],axis=1)
            self.Sconf_df = pd.DataFrame(Sconf,columns=["Sconf"])

        if self.mag == True:
            #calculate Magnetic Entropy
            if self.DLM == False:
                #Smag = np.sum(np.log(np.abs(cm_df.Mag)+1.)*cm_df.Conc,axis=1)
                Smag = kb*np.sum(np.log(np.abs(self.CM.Mag[1])+1.)*self.CM.Conc[1],axis=1)
                self.Smag_df = pd.DataFrame(Smag,columns=["Smag"])
                #self.Smag_df.reset_index(inplace=True)
                #self.Smag_df = self.Smag_df[["Smag"]]
                #self.S_df.set_index(["FN"],inplace=True)
                #self.S_df.set_index(["FN","str"],inplace=True)
                #self.ES = self.EN.join(self.S_df,on=["FN","str"])
            elif self.DLM == True:
                Smag = kb*np.sum(np.log(np.abs(self.CM.Mag[1].ix[:,1::2])+1.)*2*self.CM.Conc[1].ix[:,1::2],axis=1)
                self.Smag_df = pd.DataFrame(Smag,columns=["Smag"])
                #self.Smag_df.reset_index(inplace=True)
                #self.Smag_df = self.Smag_df[["Smag"]]

        self.CM.insert(5, 'Sconf', self.Sconf_df.loc[:,'Sconf'])
        if self.mag == True:
            self.CM.insert(6, 'Smag', self.Smag_df.loc[:,'Smag'])            

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
