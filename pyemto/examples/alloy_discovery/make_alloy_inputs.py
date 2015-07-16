import pyemto
import numpy as np
import os

latpath = "../../../" # Path do bmdl, kstr and shape directories


# each system need to have same number of alloy elements
#systems = [['Fe','Al'],['Fe','Cr']]
systems = [['Fe'],['Al']]
#concentrations = [[0.5,0.5]]
concentrations = [[1.0]]

magn = "FM" # Possible NM (Non-magnetic), FM (ferromagnetic) and 
            # DLM (Disordered local moments)

initial_sws = 3.0

# Check that initialsws is correct format
if type(initial_sws) is float:
    initial_sws = [initial_sws for x in range(3)]
elif type(initial_sws) is list:
    pass
else:
    print("ERROR: Initialsws should be float or list of 3 floats")
    exit()
if not len(initial_sws) == 3:
    print("ERROR: intialsws shoubd be a float or list of 3 floats!")
    exit()
            
            
# Sanity checks
for s in systems:
    if not len(s) == len(systems[0]):
        print("Each system need to have same number of alloy elements!")
        exit()
for c in concentrations:
    if not len(c) == len(systems[0]):
        print("Each given concetrations must have same number number as elements in system!")
        exit()

# Next check magnetic states of system and initialize splits
splits = []
if magn == "FM":
    afm = "F"
    for s in systems:
        splt = []
        for atom in s:
            if atom == "Fe":
                splt.append(2.0)
            else:
                splt.append(0.5)
        splits.append(splt)
elif magn == "DLM":
    afm = "F"
    # First duplicate each atoms and concetration
    newsystems = []
    newconcs = []
    for i in range(len(systems)):
        news = []
        newc = []
        splt = []
        for j in range(len(systems[i])):
            news.append(systems[i][j])
            news.append(systems[i][j])
            if systems[i][j] == "Fe":
                splt.append( 2.0)
                splt.append(-2.0)
            else:
                splt.append( 0.5)
                splt.append(-0.5)
        splits.append(splt)
        newsystems.append(news)
    systems = newsystems
    for c in concentrations:
        newc = []
        for conc in c:
            newc.append(conc)
            newc.append(conc)
        newconcs.append(newc)
    concentrations = newconcs
elif magn == "NM":
    afm = "P"
    for s in systems:
        splt = []
        for atom in s:
            splt.append(0.0)
        splits.append(splt)
else:
    print("Wrong magnetic state is given: " + magn)
    print("Should be one of NM, FM or DLM!")
    exit()


#print(systems)
#print(concentrations)
#print(splits)

    
#We are ready to make inputs

#for s in systems:
for si in range(len(systems)):
    s = systems[si]
    split = splits[si]
    # Create main directory 
    sname = ""
    if magn == "DLM":
        nlist = [s[i] for i in range(0,len(s),2)]            
    else:
        nlist = s
    for atom in nlist:
        sname = sname + atom
    
    # 
    # Make directories
    if not os.path.lexists(sname):
        os.makedirs(sname)

    for c in concentrations:
        # Make subdirectory for concentration
        cname = ""
        count = 0
        if magn == "DLM":
            clist = [c[i] for i in range(0,len(c),2)]
        else:
            clist = c

        for conc in clist:
            count += 1
            cname = cname +str(int(conc*1000)).zfill(4)
            
            if not count == len(clist):
                cname = cname+"-"
        apath = os.path.join(sname,cname)                
        if not os.path.lexists(apath):
            os.makedirs(apath)
        # Make subdirectory for magnetic state
        apath = os.path.join(apath,magn)
        if not os.path.lexists(apath):
            os.makedirs(apath)
        
        # Construct base jobname
        jobname = ""
        
        for i in range(len(nlist)):
            if jobname == "":
                pass
            else:
                jobname = jobname + "_"
            jobname = jobname + nlist[i].lower() + "%4.2f" % (clist[i])

        # BCC first
        alloy = pyemto.System(folder=apath)
        initialsws = initial_sws[0]
        alloy.bulk(lat='bcc',
                   jobname=jobname+"_bcc",
                   latpath=latpath,
                   sws=initialsws,
                   atoms = s,
                   concs = c,
                   splts = split,
                   afm = afm,
                   amix=0.02,
                   efmix=0.9,
                   expan='M',
                   sofc='Y',
                   xc='PBE',
                   nky=21)
        swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
        alloy.lattice_constants_batch_generate(sws=swsrange)

        # FCC second
        alloy = pyemto.System(folder=apath)
        initialsws = initial_sws[1]
        alloy.bulk(lat='fcc',
                   jobname=jobname+"_fcc",
                   latpath=latpath,
                   sws=initialsws,
                   atoms = s,
                   concs = c,
                   splts = split,
                   afm = afm,
                   amix=0.02,
                   efmix=0.9,
                   expan='M',
                   sofc='Y',
                   xc='PBE',
                   nky=21)
        swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
        alloy.lattice_constants_batch_generate(sws=swsrange)        

        # HCP last
        alloy = pyemto.System(folder=apath)
        initialsws = initial_sws[2]

        alloy.bulk(lat='hcp',
                   #jobname=jobname+"_hcp",
                   jobname=jobname, # hcp add automatically hcp string to jobname
                   latpath=latpath,
                   sws=initialsws,
                   atoms = s,
                   concs = c,
                   splts = split,
                   afm = afm,
                   amix=0.02,
                   efmix=0.9,
                   expan='M',
                   sofc='Y',
                   xc='PBE',
                   nky=21,
                   nkz=17)
        swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
        alloy.lattice_constants_batch_generate(sws=swsrange)        

