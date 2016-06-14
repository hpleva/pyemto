import matplotlib.pyplot as plt
import sys
from os import popen
import numpy as np
import pandas as pd
import pymatgen as pmg
from pyemto.utilities.utils import run_bash

#pd.options.display.max_columns = 5200
#pd.set_option('notebook_repr_html',True)
#pd.set_option('display.width', 1000)
#pd.set_option('display.max_colwidth', 100)

All = slice(None)

#Kb = 8.617332478e-5*0.073498618 #Ry/K
#Kb = 1.3806488e-23*1000#mJ/K
kb = 8.6173324E-5 #eV/K
ry2ev = 13.605698066 #eV
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
    if string == 'np.nan':
        return np.nan
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

    def __init__(self,KGRN_filenames,KFCD_filenames,suffix="prn",DLM=False):
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

        # Create a list of the filenames
        self.KGRN_filenames = run_bash('ls {}'.format(self.KGRN_filenames)).split()
        self.KFCD_filenames = run_bash('ls {}'.format(self.KFCD_filenames)).split()

        #self.xc = xc
        self.DLM = DLM
        self.EnergyColumn = [0, 1, 2, 3, 4, 5]
        self.EnergyColumnName = ["FN","SWS","ELDA","EPBE","EP07","EQNA"]
        self.StructureColumn = [0,1]
        self.StructureColumnName = ["FN","Struc"]
        self.MagneticColumn = [0, 6, 8]
        self.MagneticColumnName = ["FN","IQ","Mag"]
        self.ConcentrationColumn = [0, 3, 6, 9, 10]
        self.ConcentrationColumnName = ["FN","IQ","ITA","Conc","Elem"]
        self.DOSColumn = [0, 1]
        self.DOSColumnName = ["FN","DOSEF"]
        self.nameparser = None
        self.EN = self.Df(self.Energy(),self.EnergyColumn,self.EnergyColumnName)
        self.STR = self.Df(self.Structure(),self.StructureColumn,self.StructureColumnName)
        #self.EN = self.EN.join(self.STR.set_index(["FN"]),on=["FN"]).applymap(str2num)
        self.EN = self.STR.join(self.EN.set_index(["FN"]),on=["FN"]).applymap(str2num)

    def Energy(self):
        """
        cmd =  "grep -H TOT-{} {}".format(self.xc,self.KFCD_filenames)
        all_output = [i.split() for i in popen(cmd).readlines()]
        print(all_output)
        return all_output
        """
        #"""
        xc_list = ['LDA','PBE','P07','QNA']
        all_output = []
        
        for KFCD_file in self.KFCD_filenames:
            energy_line = ''
            for xc in xc_list:
                cmd =  "grep -H TOT-{0} {1}".format(xc,KFCD_file)
                cmd_output = popen(cmd).readlines()
                if xc == 'LDA':
                    energy_line += KFCD_file+':' + ' ' + cmd_output[0].split()[7] + ' ' + cmd_output[0].split()[4]
                # Check if QNA has not been implemented:
                elif xc == 'QNA':
                    # QNA not found:
                    if len(cmd_output) == 0:
                        energy_line += ' ' + 'np.nan'
                    # QNA was found:
                    else:
                        energy_line += ' ' + cmd_output[0].split()[4]
                else:
                    energy_line += ' ' + cmd_output[0].split()[4]
            #for i in cmd_output:
            #    all_output.append(i.split())
            all_output.append(energy_line.split())
        #print(all_output)
        return all_output
        #"""
        
    def Structure(self):
        all_output = []
        for KFCD_file in self.KFCD_filenames:
            cmd = r'grep -H mdl {} |sed "s: .*/\([^/].*\).mdl: \1:" '.format(KFCD_file)
            cmd_output = popen(cmd).readlines()
            for i in cmd_output:
                all_output.append(i.split())
        return all_output

    def Mag(self):
        all_output = []
        for KFCD_file in self.KFCD_filenames:
            cmd =  "grep -H Mag {}".format(KFCD_file) 
            cmd_output = popen(cmd).readlines()
            # Non-magnetic system return an empty list so we have to create a dummy
            if len(cmd_output) == 0:
                alt_cmd =  'grep -H \'IQ.*)\' {}'.format(KFCD_file)
                alt_cmd_output = popen(alt_cmd).readlines()
                #print(alt_cmd_output)
                cmd_output = []
                for i in alt_cmd_output:
                    #cmd_output.append(KFCD_file + ':' + '  Magnetic moment for IQ =  {0} is   0.000000 mu_B'.format(i.split()[3]))
                    cmd_output.append(KFCD_file + ':' + '  Magnetic moment for IQ =  {0} is   np.nan mu_B'.format(i.split()[3]))
            for i in cmd_output:
                all_output.append(i.split())
        #print(all_output)
        return all_output
    
        #cmd =  "grep -H Mag {}".format(KFCD_file) 
        #return [i.split() for i in popen(cmd).readlines()]

    def Concentration(self):
        all_output = []
        for KFCD_file in self.KFCD_filenames:
            cmd =  'grep -H \'IQ.*)\' {}'.format(KFCD_file) 
            cmd_output = popen(cmd).readlines()
            for i in cmd_output:
                all_output.append(i.split())
        #print(all_output)
        return all_output

        #cmd =  'grep -H \'IQ.*)\' {}'.format(self.KFCD_filenames) 
        #return [i.split() for i in popen(cmd).readlines()]

    def DOSEF(self):
	"""Extract DOS at Fermi level from the KGRN input file.
	"""
        all_output = []
	nxtsws_tag = " NXTSWS:   IQ IT ITA MMT Type  NZ ION  ELN   QTR   SPLIT  FIX  CONC"
	alat_tag = "Alat ="
	dos_tag = " Dos(Ef)    ="
	mag_tag = " Magn. mom. ="
	hop_tag = " Hopfield   ="

        for KGRN_file in self.KGRN_filenames:
	    file = open(KGRN_file,'r')
	    lines = file.readlines()
            
	    for i in range(len(lines)):
	        if nxtsws_tag in lines[i]:
		    indMin = i+2
	        elif alat_tag in lines[i]:
		    indMax = i-2
		    break
            #
            concs = np.zeros(indMax + 1 - indMin)
	    doses = np.zeros(indMax + 1 - indMin)
	    its = np.zeros(indMax + 1 - indMin)
	    #
	    ind_tmp = 0
	    for i in range(indMin,indMax+1):
	        concs[ind_tmp] = float(lines[i].split()[-1])
	        its[ind_tmp] = int(lines[i].split()[1])
	        ind_tmp += 1
	    #
	    #print("concs = ",concs)
	    #print("its = ",its)
	    #
	    num_sites = np.max(its)
	    ind_tmp = len(doses) - 1
	    # Because KGRN output file is different for non- and magnetic calculations,
	    # we have to do some additional checks to make sure we are reading the right
	    # values.
	    for i in range(len(lines)-1,indMax,-1):
	        if dos_tag in lines[i]:
		    #print(lines[i])
		    if mag_tag in lines[i+1] or hop_tag in lines[i+1]:
		        #print(lines[i+1])
		        doses[ind_tmp] = float(lines[i].split()[-1])
		        ind_tmp -= 1
		    if ind_tmp == -1:
		        break
	    #
	    #for i in range(len(doses)):
	    #    print(doses[i])
	    #
	    dos_tot = 0.0
	    for i in range(len(concs)):
	        dos_tot += concs[i]*doses[i]
	    dos_tot /= num_sites
	    #print('dos_tot = ',dos_tot)
	    dos_tot /= ry2ev
            #
            all_output.append('{0} {1}'.format(KGRN_file,dos_tot).split())
        return all_output
      
    """
    def DOSEF(self):
	# Reads total DOS, if implemented in EMTO.
        cmd =  'grep -H \' Total DOS at Fermi level =\' {}'.format(self.KGRN_filenames)
        cmd_output = popen(cmd).readlines()
        cmd_output = cmd_output[1::2]
        return [i.split() for i in cmd_output]
        
        # Component resolved DOS's
        
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
    """

    #def Formula(self):
    #    """
    #    Creates a string that represents the chemical formula of the alloy.
    #    """
        
        
    def Df(self,list,col,colname):
        return pd.DataFrame([[i[x] for x in col] for i in list],columns=colname)
         
    def create_df(self):
        self.EN.reset_index(inplace=True)
        self.EN.set_index(["FN"],inplace=True)
        #
        self.conc_df  = self.Df(self.Concentration(),self.ConcentrationColumn,self.ConcentrationColumnName)
        self.conc_df.Elem = self.conc_df.Elem.str.replace("(","")
        #
        # Extract DOS(Ef) out of KGRN output files
        self.dos_df = self.Df(self.DOSEF(),self.DOSColumn,self.DOSColumnName)

        self.mag_df = self.Df(self.Mag(),self.MagneticColumn,self.MagneticColumnName)
        self.mag_df = self.mag_df.join(self.EN,on=["FN"])
        
        if self.nameparser is None:
            pass
        else:
            fn_df = pd.DataFrame(self.mag_df.FN.str.extract(self.nameparser))
            fnc = fn_df.columns.values.tolist()
            fnc.append("FN")
            fnc.append("Struc")
            self.mag_df = self.mag_df.join(fn_df)
            
        #self.cn = cn_df
        #self.mag_df = mag_df
        self.main_df  = self.conc_df.join(self.mag_df,rsuffix='r')
        if (self.main_df.FN != self.main_df.FNr).any():
            print "big problem !!!"
            raise ValueError
        self.main_df.drop(["FNr","IQr"],axis=1,inplace=True)
        
        self.main_df = pd.pivot_table(self.main_df.applymap(str2num),index=["FN","Struc","SWS","ELDA","EPBE","EP07","EQNA"],values=["Mag","Conc"],columns=["IQ","ITA","Elem"])
        #may be a bug here, we can not get Elem any more, but it seem no harm.
        
        #if self.nameparser is None:
        #    cm_df = pd.pivot_table(cm_df.applymap(str2num),index=["FN","Struc","SWS","Etot"],values=["Mag","CC","Elem"],columns=["IQ","ITA"])
        #else:
        #    cm_df = pd.pivot_table(cm_df.applymap(str2num),index=fnc,values=["Mag","CC","Elem"],columns=["IQ","ITA"])

        # Some tricks
        self.main_df = self.main_df.reset_index()
        # Get rid of the ':' that grep puts at the end of the filename:
        self.main_df.FN = self.main_df.FN.str.replace(":","")
        
        #"""
        # Calculate configurational entropy
        if self.DLM == False:
            Sconf = -kb*np.sum(np.log(self.main_df.Conc[1])*self.main_df.Conc[1],axis=1)
            self.Sconf_df = pd.DataFrame(Sconf,columns=["Sconf"])
            #self.Sconf_df.reset_index(inplace=True)
            #self.Sconf_df = self.Sconf_df[["Sconf"]]
        elif self.DLM == True:
            Sconf = -kb*np.sum(np.log(2*self.main_df.Conc[1].ix[:,1::2])*2*self.main_df.Conc[1].ix[:,1::2],axis=1)
            self.Sconf_df = pd.DataFrame(Sconf,columns=["Sconf"])
        #"""
        #"""
        #calculate Magnetic Entropy
        if self.DLM == False:
            #Smag = np.sum(np.log(np.abs(cm_df.Mag)+1.)*cm_df.Conc,axis=1)
            Smag = kb*np.sum(np.log(np.abs(self.main_df.Mag[1])+1.)*self.main_df.Conc[1],axis=1)
            self.Smag_df = pd.DataFrame(Smag,columns=["Smag"])
            #self.Smag_df.reset_index(inplace=True)
            #self.Smag_df = self.Smag_df[["Smag"]]
            #self.S_df.set_index(["FN"],inplace=True)
            #self.S_df.set_index(["FN","Struc"],inplace=True)
            #self.ES = self.EN.join(self.S_df,on=["FN","Struc"])
        elif self.DLM == True:
            Smag = kb*np.sum(np.log(np.abs(self.main_df.Mag[1].ix[:,1::2])+1.)*2*self.main_df.Conc[1].ix[:,1::2],axis=1)
            self.Smag_df = pd.DataFrame(Smag,columns=["Smag"])
            #self.Smag_df.reset_index(inplace=True)
            #self.Smag_df = self.Smag_df[["Smag"]]
        #"""

        self.main_df.insert(7, 'Sconf', self.Sconf_df.loc[:,'Sconf'])
        self.main_df.insert(8, 'Smag', self.Smag_df.loc[:,'Smag'])
        self.main_df.insert(9, 'DOSEF', self.dos_df.loc[:,'DOSEF'])

        # Construct the chemical formula of the system:
        #self.formula_df
            
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
