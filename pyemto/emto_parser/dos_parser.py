import numpy as np

class Atom(object):
    def __init__(self, i, label, arr1, arr2):
        self.number = i
        self.label = label
        self.e = arr1[:,0]
        self.dos = arr1[:,1]
        self.sdos = arr1[:,2]
        self.pdos = arr1[:,3]
        self.ddos = arr1[:,4]
        self.fdos = arr1[:,5]
        #
        self.nos = arr2[:,1]
        self.snos = arr2[:,2]
        self.pnos = arr2[:,3]
        self.dnos = arr2[:,4]
        self.fnos = arr2[:,5]
        
class Element(object):
    def __init__(self, dos, nos, sdos, pdos, ddos, fdos, snos, pnos, dnos, fnos):
        self.dos = dos
        self.nos = nos
        self.sdos = sdos
        self.pdos = pdos
        self.ddos = ddos
        self.fdos = fdos
        self.snos = snos
        self.pnos = pnos
        self.dnos = dnos
        self.fnos = fnos
        
        
class Get_DOS(object):
    """Parser for EMTO DOS files"""
    def __init__(self, fn):
        self.fn = fn
        self.data = None
        with open(self.fn, 'r') as self.data:
            self.lines = [line.rstrip() for line in self.data]
        self.e = None
        # Total DOS and NOS
        self.dos = None
        self.nos = None
        # Lattice data
        self.nq = 0
        self.labels = []
        # Collect info and find breakpoints
        last_line = ''
        self.startpoints = []
        self.endpoints = []
        for i, line in enumerate(self.lines):
            if 'Sublattice' in line:
                tmp = line.split()[3]
                self.nq += 1
                self.labels.append(tmp)
            try:
                tmp1 = float(line.split()[0])
                try:
                    tmp2 = float(last_line.split()[0])
                except:
                    self.startpoints.append(i)
                    #print(last_line)
                    #print(line)
                    #print()
            except:
                try:
                    tmp2 = float(last_line.split()[0])
                    self.endpoints.append(i-1)
                except:
                    pass
            last_line = line
        # Get list of unique elements
        self.unique = []
        self.unique_ind = []
        for label in self.labels:
            if label not in self.unique:
                self.unique.append(label)
                self.unique_ind.append(np.where(np.array(self.labels) == label)[0])
        # Get actual data
        self.atoms = []
        self.tmp = np.array([l.split() for l in self.lines[self.startpoints[0]:self.endpoints[0]+1]], dtype=np.float)
        self.e = self.tmp[:,0]
        self.dos = self.tmp[:,1]
        self.nos = self.tmp[:,2]
        for i in range(self.nq):
            ind1 = (i+1)*2 - 1
            ind2 = (i+1)*2
            tmp1 = np.array([l.split() for l in self.lines[self.startpoints[ind1]:self.endpoints[ind1]+1]], dtype=np.float)
            tmp2 = np.array([l.split() for l in self.lines[self.startpoints[ind2]:self.endpoints[ind2]+1]], dtype=np.float)
            self.atoms.append(Atom(i, self.labels[i], tmp1, tmp2))
        # Get total s,p,d,f DOS and NOS
        self.sdos = np.zeros_like(self.e)
        self.pdos = np.zeros_like(self.e)
        self.ddos = np.zeros_like(self.e)
        self.fdos = np.zeros_like(self.e)
        self.snos = np.zeros_like(self.e)
        self.pnos = np.zeros_like(self.e)
        self.dnos = np.zeros_like(self.e)
        self.fnos = np.zeros_like(self.e)
        for atom in self.atoms:
            self.sdos += atom.sdos
            self.pdos += atom.pdos
            self.ddos += atom.ddos
            self.fdos += atom.fdos
            self.snos += atom.snos
            self.pnos += atom.pnos
            self.dnos += atom.dnos
            self.fnos += atom.fnos
        # Get element-wise total DOS and NOS
        self.elems = {}
        for i, u in enumerate(self.unique):
            dos = np.zeros_like(self.e)
            nos = np.zeros_like(self.e)
            sdos = np.zeros_like(self.e)
            pdos = np.zeros_like(self.e)
            ddos = np.zeros_like(self.e)
            fdos = np.zeros_like(self.e)
            snos = np.zeros_like(self.e)
            pnos = np.zeros_like(self.e)
            dnos = np.zeros_like(self.e)
            fnos = np.zeros_like(self.e)
            for j in self.unique_ind[i]:
                dos += self.atoms[j].dos
                nos += self.atoms[j].nos
                sdos += self.atoms[j].sdos
                pdos += self.atoms[j].pdos
                ddos += self.atoms[j].ddos
                fdos += self.atoms[j].fdos
                snos += self.atoms[j].snos
                pnos += self.atoms[j].pnos
                dnos += self.atoms[j].dnos
                fnos += self.atoms[j].fnos
            self.elems[self.unique[i]] = Element(dos, nos, sdos, pdos, ddos, fdos, snos, pnos, dnos, fnos)
        
    def print_file(self):
        for i, line in enumerate(self.lines):
            #if '#' in line:
                print(i, line)
