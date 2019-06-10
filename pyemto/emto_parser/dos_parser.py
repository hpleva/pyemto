import numpy as np

class Atom(object):
    def __init__(self, i, label, spin, arr1, arr2):
        self.sublattice = i
        self.label = label
        self.spin = spin
        self.e = arr1[:,0]
        if self.spin == 'up':
            self.dos = arr1[:,1]
            self.sdos = arr1[:,2]
            self.pdos = arr1[:,3]
            self.ddos = arr1[:,4]
            self.fdos = arr1[:,5]
        else:
            self.dos = -np.abs(arr1[:,1])
            self.sdos = -np.abs(arr1[:,2])
            self.pdos = -np.abs(arr1[:,3])
            self.ddos = -np.abs(arr1[:,4])
            self.fdos = -np.abs(arr1[:,5])
        #
        self.nos = arr2[:,1]
        self.snos = arr2[:,2]
        self.pnos = arr2[:,3]
        self.dnos = arr2[:,4]
        self.fnos = arr2[:,5]

class Element(object):
    def __init__(self, e, dos, nos, sdos, pdos, ddos, fdos, snos, pnos, dnos, fnos):
        self.e = e
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
    def __init__(self, fn):
        self.fn = fn
        self._data = None
        with open(self.fn, 'r') as self._data:
            self._lines = [line.rstrip() for line in self._data]

        _EF_dim = 100
        self._EFx = np.zeros(_EF_dim)
        self._EFy = np.linspace(-100, 100, _EF_dim)

        # First get a list of all breakpoints
        last_line = ''
        self.startpoints = []
        self.endpoints = []
        for i, line in enumerate(self._lines):
            try:
                tmp1 = float(line.split()[0])
                try:
                    tmp2 = float(last_line.split()[0])
                except:
                    self.startpoints.append(i)
            except:
                try:
                    tmp2 = float(last_line.split()[0])
                    self.endpoints.append(i-1)
                except:
                    pass
            last_line = line

        bp_index = 0
        self.mag = False
        self.labels = []
        spin_situation_is_clear = False
        self._structure = {'total':{}, 'sublat':{}}
        for i, line in enumerate(self._lines):
            if 'Total DOS' in line:
                if 'DOSUP+DOWN' in line or 'DOSUP' in line:
                    spin = 'up'
                elif 'DOSDOWN' in line:
                    spin = 'down'
                    if not spin_situation_is_clear:
                        self.mag = True
                        spin_situation_is_clear = True
                self._structure['total'][spin] = {'start': self.startpoints[bp_index],
                                                  'end': self.endpoints[bp_index]}
                bp_index += 1
                
            elif  'Sublattice' in line:
                sublat_ind = int(line.split()[1])
                sublat_elem = line.split()[3]
                sublat_spin = line.split()[5].lower()
                self.labels.append(sublat_elem)
                if sublat_spin == 'up+down':
                    sublat_spin = 'up'
                if sublat_spin not in self._structure['sublat']:
                    self._structure['sublat'][sublat_spin] = {}
                if sublat_ind not in self._structure['sublat'][sublat_spin]:
                    self._structure['sublat'][sublat_spin][sublat_ind] = {}
                self._structure['sublat'][sublat_spin][sublat_ind][sublat_elem] = \
                        {'start': self.startpoints[bp_index],
                         'end': self.endpoints[bp_index],
                         'tnos_start': self.startpoints[bp_index+1],
                         'tnos_end': self.endpoints[bp_index+1]}
                bp_index += 2

        # Get list of unique elements
        self.unique = []
        self.unique_ind = []
        for label in self.labels:
            if label not in self.unique:
                self.unique.append(label)
                self.unique_ind.append(np.where(np.array(self.labels) == label)[0])
                
        # Now that we have the structure of the data figured out,
        # we can start extracting the data.
        
        # Total values
        start = self._structure['total']['up']['start']
        end = self._structure['total']['up']['end']
        tmp = np.array([l.split() for l in self._lines[start:end+1]], dtype=np.float)
        dim = (len(tmp[:,0]), 2)
        self.e = np.zeros(len(tmp[:,0]))
        self.dos = np.zeros(dim)
        self.nos = np.zeros(dim)
        self.e = tmp[:,0]
        self.dos[:,0] = tmp[:,1]
        self.nos[:,0] = tmp[:,2]
        if self.mag:
            start = self._structure['total']['down']['start']
            end = self._structure['total']['down']['end']
            tmp = np.array([l.split() for l in self._lines[start:end+1]], dtype=np.float)
            self.dos[:,1] = -np.abs(tmp[:,1])
            self.nos[:,1] = tmp[:,2]

        self.atoms = []
        for spin in self._structure['sublat'].keys():
            for key in self._structure['sublat'][spin].keys():
                for key2 in self._structure['sublat'][spin][key]:
                    dos_start = self._structure['sublat'][spin][key][key2]['start']
                    dos_end = self._structure['sublat'][spin][key][key2]['end']
                    nos_start = self._structure['sublat'][spin][key][key2]['tnos_start']
                    nos_end = self._structure['sublat'][spin][key][key2]['tnos_end']
                    tmp1 = np.array([l.split() for l in self._lines[dos_start:dos_end+1]], dtype=np.float)
                    tmp2 = np.array([l.split() for l in self._lines[nos_start:nos_end+1]], dtype=np.float)
                    # Concentrations per sublattice element have already been taken care of in KGRN.
                    self.atoms.append(Atom(key, key2, spin, tmp1, tmp2))

        # Get total s,p,d,f DOS and NOS
        self.sdos = np.zeros_like(self.dos)
        self.pdos = np.zeros_like(self.dos)
        self.ddos = np.zeros_like(self.dos)
        self.fdos = np.zeros_like(self.dos)
        self.snos = np.zeros_like(self.dos)
        self.pnos = np.zeros_like(self.dos)
        self.dnos = np.zeros_like(self.dos)
        self.fnos = np.zeros_like(self.dos)
        for atom in self.atoms:
            if atom.spin == 'up':
                ind = 0
            else:
                ind = 1
            self.sdos[:,ind] += atom.sdos
            self.pdos[:,ind] += atom.pdos
            self.ddos[:,ind] += atom.ddos
            self.fdos[:,ind] += atom.fdos
            self.snos[:,ind] += atom.snos
            self.pnos[:,ind] += atom.pnos
            self.dnos[:,ind] += atom.dnos
            self.fnos[:,ind] += atom.fnos
            
        # Get element-wise total DOS and NOS
        self.elems = {}
        for i, u in enumerate(self.unique):
            dos = np.zeros_like(self.dos)
            nos = np.zeros_like(self.dos)
            sdos = np.zeros_like(self.dos)
            pdos = np.zeros_like(self.dos)
            ddos = np.zeros_like(self.dos)
            fdos = np.zeros_like(self.dos)
            snos = np.zeros_like(self.dos)
            pnos = np.zeros_like(self.dos)
            dnos = np.zeros_like(self.dos)
            fnos = np.zeros_like(self.dos)
            for j in self.unique_ind[i]:
                if self.atoms[j].spin == 'up':
                    ind = 0
                else:
                    ind = 1
                dos[:,ind] += self.atoms[j].dos
                nos[:,ind] += self.atoms[j].nos
                sdos[:,ind] += self.atoms[j].sdos
                pdos[:,ind] += self.atoms[j].pdos
                ddos[:,ind] += self.atoms[j].ddos
                fdos[:,ind] += self.atoms[j].fdos
                snos[:,ind] += self.atoms[j].snos
                pnos[:,ind] += self.atoms[j].pnos
                dnos[:,ind] += self.atoms[j].dnos
                fnos[:,ind] += self.atoms[j].fnos
            self.elems[self.unique[i]] = Element(self.e, dos, nos, sdos, pdos, ddos, fdos, snos, pnos, dnos, fnos)

    def print_file(self):
        for i, line in enumerate(self._lines):
            print(i, line)

    @property
    def EFx(self):
        return self._EFx

    @property
    def EFy(self):
        return self._EFy
