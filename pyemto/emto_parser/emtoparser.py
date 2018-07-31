import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import six
from pyemto.utilities.utils import run_bash

try:
    import pandas as pd
except ImportError:
    raise ImportError('EMTOPARSER requires pandas>=0.20.3 to be installed!')

# Set up some constants
ry2ev = 13.605698066 #eV
kb = 8.6173324E-5/ry2ev #Ry/K
bohr2angstrom = 0.529177249

# Loads periodic table data from a json file
with open(os.path.join(os.path.dirname(__file__), "periodic_table.json"), "rt") as f:
    periodic_table = pd.read_json(f)

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
    if not isinstance(string, str):
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

            
class EMTOPARSER(object):
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

        # In python3 run_bash returns bytestrings, so first we make sure they are decoded
        # to regular strings.
        for i in range(len(self.KGRN_filenames)):
            self.KGRN_filenames[i] = self.KGRN_filenames[i].decode()
        for i in range(len(self.KFCD_filenames)):
            self.KFCD_filenames[i] = self.KFCD_filenames[i].decode()

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
        self.NQColumn = [0, 3]
        self.NQColumnName = ["FN","NQ"]
        self.FormulaColumn = [0,1]
        self.FormulaColumnName = ["Formula","Mav"]
        self.StatusColumn = [0,1]
        self.StatusColumnName = ["FN","Status"]
        self.COAColumn = [0,9]
        self.COAColumnName = ["FN","COA"]
        #self.MavColumn = [0]
        #self.MavColumnName = ["Mav"]
        #
        self.nameparser = None
        self.EN = self.Df(self.Energy(),self.EnergyColumn,self.EnergyColumnName).applymap(str2num)
        self.STR = self.Df(self.Structure(),self.StructureColumn,self.StructureColumnName)
        self.STATUS = self.Df(self.get_Status(),self.StatusColumn,self.StatusColumnName)
        self.COA = self.Df(self.get_COA(),self.COAColumn,self.COAColumnName)
        #self.EN = self.EN.join(self.STR.set_index(["FN"]),on=["FN"]).applymap(str2num)
        self.EN = self.STR.join(self.EN.set_index(["FN"]), on=["FN"]).applymap(str2num)
        self.EN.insert(1, 'Status',   self.STATUS.loc[:,'Status'])
        self.EN.insert(3, 'COA', self.COA.loc[:,'COA'])

    def Energy(self):
        """
        Extract Total energies from the KFCD output file
        """
        xc_list = ['LDA','PBE','P07','QNA']
        all_output = []

        for KFCD_file in self.KFCD_filenames:
            energy_line = ''
            for xc in xc_list:
                cmd =  "grep -H TOT-{0} {1}".format(xc, KFCD_file)
                cmd_output = os.popen(cmd).readlines()
                if xc == 'LDA':
                    #try:
                    #    energy_line += KFCD_file+':' + ' ' + cmd_output[0].split()[7] + ' ' + cmd_output[0].split()[4]
                    #except:
                    #    print(KFCD_file)
                    #    quit()
                    # If KGRN didnt finish properly cmd_output is empty, causing a crash right here!
                    if len(cmd_output) == 0:
                        print(KFCD_file + ' didn\'t finish properly!')
                        energy_line += KFCD_file+':' + ' ' + '0.0' + ' ' + '0.0'
                    else:
                        energy_line += KFCD_file+':' + ' ' + cmd_output[0].split()[7] + ' ' + cmd_output[0].split()[4]
                # Check if QNA has not been implemented:
                elif xc == 'QNA':
                    # QNA not found:
                    if len(cmd_output) == 0:
                        #energy_line += ' ' + 'np.nan'
                        energy_line += ' ' + '0.0'
                    # QNA was found:
                    else:
                        energy_line += ' ' + cmd_output[0].split()[4]
                else:
                    if len(cmd_output) == 0:
                        energy_line += ' ' + '0.0'
                    else:
                        energy_line += ' ' + cmd_output[0].split()[4]
            #for i in cmd_output:
            #    all_output.append(i.split())
            all_output.append(energy_line.split())
        #for ii in all_output:
        #    print(ii)
        return all_output
        #"""

    def Structure(self):
        """ Extract structure name from the BDML file.
        """
        all_output = []
        for KFCD_file in self.KFCD_filenames:
            cmd = r'grep -H mdl {} |sed "s: .*/\([^/].*\).mdl: \1:" '.format(KFCD_file)
            cmd_output = os.popen(cmd).readlines()
            for i in cmd_output:
                all_output.append(i.split())
        return all_output

    def Mag(self):
        """  Extract magnetic momets from KFCD output 
        """
        all_output = []
        for KFCD_file in self.KFCD_filenames:
            cmd =  "grep -H Mag {}".format(KFCD_file)
            cmd_output = os.popen(cmd).readlines()
            # Non-magnetic systems return an empty list so we have to create a dummy:
            if len(cmd_output) == 0:
                alt_cmd =  'grep -H \'IQ.*)\' {}'.format(KFCD_file)
                alt_cmd_output = os.popen(alt_cmd).readlines()
                # Crashed calculations will also return an empty list even here:
                if len(alt_cmd_output) == 0:
                    cmd_output = [KFCD_file+':           Magnetic moment for IQ =  1 is   0.000000 mu_B\n']
                else:
                    cmd_output = []
                    for i in alt_cmd_output:
                        cmd_output.append(KFCD_file + ':' + '  Magnetic moment for IQ =  {0} is   0.000000 mu_B'.format(i.split()[3]))
                        #cmd_output.append(KFCD_file + ':' + '  Magnetic moment for IQ =  {0} is   np.nan mu_B'.format(i.split()[3]))
            #print(cmd_output)
            for i in cmd_output:
                all_output.append(i.split())
        #print(all_output)
        return all_output

        #cmd =  "grep -H Mag {}".format(KFCD_file)
        #return [i.split() for i in os.popen(cmd).readlines()]

    def Concentration(self):
        """ Extract elemental concentration for each site from KFCD output
        """
        all_output = []
        for KFCD_file in self.KFCD_filenames:
            cmd =  'grep -H \'IQ.*)\' {}'.format(KFCD_file)
            cmd_output = os.popen(cmd).readlines()
            # Sanity check for crashed calculations:
            if len(cmd_output) == 0:
                cmd_output = [KFCD_file+':    IQ =   1 ITA =  1 CONC =  0.000000  (Er  )\n']
            #print(cmd_output)
            for i in cmd_output:
                all_output.append(i.split())
        #print(all_output)
        return all_output

        #cmd =  'grep -H \'IQ.*)\' {}'.format(self.KFCD_filenames)
        #return [i.split() for i in os.popen(cmd).readlines()]

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
            for i in range(indMin, indMax+1):
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
            dos_tot = 0.0
            for i in range(len(concs)):
                dos_tot += concs[i]*doses[i]
            dos_tot /= num_sites
            #
            all_output.append('{0} {1}'.format(KGRN_file,dos_tot).split())
        return all_output

    def get_Formula_and_Mav(self):
        """
        Creates a string that represents the chemical formula of the alloy.
        """
        all_output = []
        index_max = self.main_df.shape[0]
        #
        for index in range(index_max):
            elems = self.main_df[index:index+1].Elem.values.flatten()
            concs = self.main_df[index:index+1].Conc.values.flatten()
            nq_tmp = self.main_df.NQ[index]
            elem_len = len(elems)
            formula_tmp = {}
            #print(nq_tmp, elems)
            for i in range(elem_len):
                if elems[i] == None:
                    pass
                else:
                    if elems[i] in formula_tmp:
                        formula_tmp[elems[i]] += concs[i]/nq_tmp
                    else:
                        formula_tmp[elems[i]] = concs[i]/nq_tmp
            string_tmp = ''
            mav_tmp = 0.0
            # iteritems method return stuff in the order it is
            # organized internally inside the dictionary.
            # To get alphabetical order we can use sorted() method:
            for key, value in sorted(six.iteritems(formula_tmp)):
                #print(key,value)
                string_tmp += '{0:2}{1:5.3f}'.format(key,value)
                mav_tmp += periodic_table[key]['Atomic mass']*value
            #print(string_tmp)
            all_output.append([string_tmp,mav_tmp])
        return all_output

    def get_NQ(self):
        """ 
        Extract number of atom positions from KFCD
        """
        all_output = []
        for KFCD_file in self.KFCD_filenames:
            cmd = 'grep -H NQ {}'.format(KFCD_file)
            cmd_output = os.popen(cmd).readlines()
            # Sanity check for crashed calculations:
            if len(cmd_output) == 0:
                cmd_output = [KFCD_file+":           NQ  =  1 NL  =  0 NLM = 81 NS  =  0 CPA_ALPHA =  0.000000\n"]
            for i in cmd_output:
                all_output.append(i.split())
        return all_output

    def get_COA(self):
        """
        Extract C over A ration from the KFCD output
        """
        
        all_output = []
        for KFCD_file in self.KFCD_filenames:
            cmd = 'grep -H \'BSZ( 3)\' {}'.format(KFCD_file)
            cmd_output = os.popen(cmd).readlines()
            # Sanity check for crashed calculations:
            if len(cmd_output) == 0:
                cmd_output = [KFCD_file+':           BSX( 3)=   0.00000 BSY( 3)=   0.00000 BSZ( 3)=   0.00000\n']
            for i in cmd_output:
                all_output.append(i.split())
        return all_output

    def get_Status(self):
        """
        Extract status of calculation from KGRN output
        """
        all_output = []
        for KGRN_file in self.KGRN_filenames:
            cmd = 'grep -H \' Converged in\' {}'.format(KGRN_file)
            cmd_output = os.popen(cmd).readlines()
            # Check if not converged:
            if len(cmd_output) == 0:
                cmd = 'grep -H \' Not converged\' {}'.format(KGRN_file)
                cmd_output = os.popen(cmd).readlines()
                # Check if crashed somehow:
                if len(cmd_output) == 0:
                    cmd = 'grep -H \'Fermi level not found.\' {}'.format(KGRN_file)
                    cmd_output = os.popen(cmd).readlines()
                    if len(cmd_output) != 0:
                        cmd_output = [cmd_output[-1].rstrip('\n').replace(' ', '_').replace(':_', ': ')]
                    else:
                        cmd_output = [KGRN_file+':' + ' ' + 'Crashed_somehow!!!']
                else:
                    cmd_output = [KGRN_file+':' + ' ' + 'Not_converged!!!']
            else:
                split_tmp = cmd_output[0].split()
                if len(split_tmp) == 7:
                    cmd_output = [split_tmp[0] + ' ' +\
                                  split_tmp[1] + '_' +\
                                  split_tmp[2] + '_' +\
                                  split_tmp[3] + '_' +\
                                  split_tmp[4] + '_' +\
                                  split_tmp[5] + '_' +\
                                  split_tmp[6]]
                elif len(split_tmp) == 8:
                    cmd_output = [split_tmp[0] + ' ' +\
                                  split_tmp[1] + '_' +\
                                  split_tmp[2] + '_' +\
                                  split_tmp[3] + '_' +\
                                  split_tmp[4] + '_' +\
                                  split_tmp[5] + '_' +\
                                  split_tmp[6] + '_' +\
                                  split_tmp[7]]
            all_output.append(cmd_output[0].split())
        return all_output

    def Df(self,list,col,colname):
        return pd.DataFrame([[i[x] for x in col] for i in list],columns=colname)

    def create_df(self):
        """
        Create pandas DataFrame from extracted data
        """
        self.EN.reset_index(inplace=True)
        self.EN.set_index(["FN"], inplace=True)
        #
        self.conc_df  = self.Df(self.Concentration(), self.ConcentrationColumn, self.ConcentrationColumnName)
        self.conc_df.Elem = self.conc_df.Elem.str.replace("(", "")
        #
        # Extract DOS(Ef) out of KGRN output files
        self.dos_df = self.Df(self.DOSEF(), self.DOSColumn, self.DOSColumnName).applymap(str2num)

        # Extract magnetic moments
        self.mag_df = self.Df(self.Mag(), self.MagneticColumn, self.MagneticColumnName)
        self.mag_df = self.mag_df.join(self.EN, on=["FN"])

        if self.nameparser is None:
            pass
        else:
            fn_df = pd.DataFrame(self.mag_df.FN.str.extract(self.nameparser))
            fnc = fn_df.columns.values.tolist()
            fnc.append("FN")
            fnc.append("Struc")
            self.mag_df = self.mag_df.join(fn_df)

        self.main_df  = self.conc_df.join(self.mag_df, rsuffix='r')
        if (self.main_df.FN != self.main_df.FNr).any():
            print("big problem !!!")
            raise ValueError
        self.main_df.drop(["FNr", "IQr"], axis=1, inplace=True)

        self.main_df = pd.pivot_table(self.main_df.applymap(str2num),
        index=["index", "FN", "Status", "Struc", "COA", "SWS", "ELDA", "EPBE", "EP07", "EQNA"],
        values=["Mag", "Conc", "Elem"], columns=["IQ", "ITA"], aggfunc = lambda x: max(x)
        )

        self.main_df = self.main_df.applymap(str2num)
        self.main_df = self.main_df.reset_index()
        self.main_df.drop(["index"], axis=1, inplace=True)
        # Get rid of the ':' that grep puts at the end of the filename:
        self.main_df.FN = self.main_df.FN.str.replace(":", "")

        self.nq_df = self.Df(self.get_NQ(), self.NQColumn, self.NQColumnName).applymap(str2num)

        # Compute volume per atom
        vol_per_atom = 4.0/3*np.pi*(self.main_df.SWS*bohr2angstrom)**3
        self.main_df.insert(4, 'VOL', vol_per_atom)

        # Calculate configurational entropy
        if self.DLM == False:
            Sconf = -kb*np.sum(np.log(self.main_df.Conc[1])*self.main_df.Conc[1],axis=1)
            self.Sconf_df = pd.DataFrame(Sconf,columns=["Sconf"])
        elif self.DLM == True:
            Sconf = -kb*np.sum(np.log(2*self.main_df.Conc[1].ix[:,::2])*2*self.main_df.Conc[1].ix[:,::2],axis=1)
            self.Sconf_df = pd.DataFrame(Sconf,columns=["Sconf"])
        #calculate Magnetic Entropy
        if self.DLM == False:
            Smag = kb*np.sum(np.log(np.abs(self.main_df.Mag[1])+1.)*self.main_df.Conc[1],axis=1)
            self.Smag_df = pd.DataFrame(Smag,columns=["Smag"])
        elif self.DLM == True:
            Smag = kb*np.sum(np.log(np.abs(self.main_df.Mag[1].ix[:,::2])+1.)*2*self.main_df.Conc[1].ix[:,::2],axis=1)
            self.Smag_df = pd.DataFrame(Smag,columns=["Smag"])

        insert_index = 10
        self.main_df.insert(insert_index,   'Sconf',   self.Sconf_df.loc[:,'Sconf'])
        self.main_df.insert(insert_index+1, 'Smag',    self.Smag_df.loc[:,'Smag'])
        self.main_df.insert(insert_index+2, 'DOSEF',   self.dos_df.loc[:,'DOSEF'])
        self.main_df.insert(insert_index+3, 'NQ',      self.nq_df.loc[:,'NQ'])

        # Construct the chemical formula and the average atomic mass of the system:
        self.formula_df = self.Df(self.get_Formula_and_Mav(),self.FormulaColumn,self.FormulaColumnName).applymap(str2num)
        self.main_df.insert(insert_index+3, 'Formula', self.formula_df.loc[:,'Formula'])
        self.main_df.insert(insert_index+4, 'Mav', self.formula_df.loc[:,'Mav'])

        # Calculate average atomic mass (which is needed when Debye temperature is calculated.)
        #self.mav_df = self.Df(self.get_Mav(),self.MavColumn,self.MavColumnName)


if __name__=='__main__' :
    test = EMTOPARSER("kfcd_output/al0.20_co0.20_cr0.20_fe0.20_ni0.20_bcc_2.550000",suffix='prn')
    print(test.Energy())
    print(test.Structure())
    print(test.Concentration())
    test.MagTables()
    print(test.CM)
