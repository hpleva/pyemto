import numpy as np
import sys

class emto_dos:
    def __init__(self,dos_path):
        """Takes an EMTO DOS output file and parses the data into numpy arrays for easy plotting etc.
        """

        self.dos_path = dos_path
        self.dos_lines = None
        self.dos_auth_str = '     Total DOS and NOS and partial (IT) DOSUP+DOWN'

        self.E_total_str = '     E     Total     1'
        self.E_total_sublattice_str = '     E     Total       s       p       d       f'
        self.tnos_str = '  TNOS'
        
        self.total_dos_nos_array = None
        self.total_dos_nos_index_beg = None
        self.total_dos_nos_index_end = None

        self.sublattice_dos_spdf_array = []
        self.sublattice_dos_spdf_index_beg = []
        self.sublattice_dos_spdf_index_end = []

        self.sublattice_tnos_spdf_array = []
        self.sublattice_tnos_spdf_index_beg = []
        self.sublattice_tnos_spdf_index_end = []

        self.sublattice_info = []
        
        with open(dos_path, "rt") as dos_file:
            self.dos_lines = [line.rstrip() for line in dos_file]
            #print(self.dos_lines)
            #for line in self.dos_lines:
            #    print(line)
            
            # Confirm to some extent that this is a valid EMTO DOS output file:
            if self.dos_auth_str not in self.dos_lines:
                raise Exception('File {0} does not look like a valid EMTO DOS file!\n'.format(dos_path)+
                                '           File does not contain the usual header \'{0}\'!'.format(self.dos_auth_str))
            
            # Extract line indexes of the beginnings and ends of data blocks:
            index = 0
            first_sublattice = True
            for line in self.dos_lines:
                if self.E_total_str in line:
                    self.total_dos_nos_index_beg = index+2
                elif self.E_total_sublattice_str in line:
                    if first_sublattice == True:
                        self.total_dos_nos_index_end = index-4
                        first_sublattice = False
                    else:
                        self.sublattice_tnos_spdf_index_end.append(index-5)
                    self.sublattice_dos_spdf_index_beg.append(index+2)
                    sublattice_line = self.dos_lines[index-2]
                    sublattice_line_list = sublattice_line.split()
                    #print(sublattice_line)
                    self.sublattice_info.append([sublattice_line_list[1],sublattice_line_list[3]])
                elif self.tnos_str in line:
                    self.sublattice_tnos_spdf_index_beg.append(index+2)
                    self.sublattice_dos_spdf_index_end.append(index-2)
                index += 1
            # Get the last TNOS index
            for i in xrange(index-1,0,-1):
                if self.dos_lines[i] != '':
                    self.sublattice_tnos_spdf_index_end.append(i)

            # Now populate the arrays
            # E, total DOS, total NOS, THEN total DOS per sublattice (spin up + spin down)
            self.total_dos_nos_array = np.loadtxt(self.dos_lines[self.total_dos_nos_index_beg:self.total_dos_nos_index_end+1])

            #"""
            for i in range(len(self.sublattice_info)):
                #      E     Total DOS       s       p       d       f
                self.sublattice_dos_spdf_array.append(np.loadtxt(self.dos_lines[self.sublattice_dos_spdf_index_beg[i]:self.sublattice_dos_spdf_index_end[i]+1]))
                #      E     Total NOS       s       p       d       f
                self.sublattice_tnos_spdf_array.append(np.loadtxt(self.dos_lines[self.sublattice_tnos_spdf_index_beg[i]:self.sublattice_tnos_spdf_index_end[i]+1]))
            #"""
                
            #for row in self.sublattice_tnos_spdf_array[3]:
            #for row in self.dos_lines[self.sublattice_dos_spdf_index_beg[i]:self.sublattice_dos_spdf_index_end[i]+1]:
            #    print(row)
            #print(test)
            #i = 3
            #print(self.dos_lines[self.sublattice_dos_spdf_index_beg[i]:self.sublattice_dos_spdf_index_end[i]+1])
            #print(self.dos_lines[self.sublattice_tnos_spdf_index_beg[i]:self.sublattice_tnos_spdf_index_end[i]+1])

            #self.total_dos_nos_array = np.loadtxt(,)
                    
            
    def plot_dos(self):
        import pylab
        pylab.plot([0.0,0.0],[-10,100],'--',linewidth=2,color='red')
        #total dos
        #pylab.plot(self.total_dos_nos_array[:,0],self.total_dos_nos_array[:,1])
        #component dos
        pylab.plot(self.sublattice_dos_spdf_array[0][:,0],self.sublattice_dos_spdf_array[0][:,1])
        pylab.plot(self.sublattice_dos_spdf_array[1][:,0],self.sublattice_dos_spdf_array[1][:,1],'--')
        pylab.xlim(-0.6,0.05)
        pylab.ylim(-1.0,40)
        pylab.title(self.dos_path)
        pylab.show()
