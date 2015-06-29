import pyemto
import numpy as np
import os

latpath = "../../../../" # Path do bmdl, kstr and shape directories


# each system need to have same number of alloy elements
#systems = [['Fe','Al'],['Fe','Cr']]
#systems = [['Fe'],['Al']]
systems = [['Al']]
#concentrations = [[0.5,0.5]]
concentrations = [[1.0]]

magn = "NM" # Possible NM (Non-magnetic), FM (ferromagnetic) and 
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

results = []
#We are ready to make inputs
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
        sc_res = []
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
        finalname = jobname + "_final"

        # BCC first
        alloy = pyemto.System(folder=apath)
        initialsws = initial_sws[0] # We need some clever way to get this
        alloy.bulk(lat='bcc', jobname=jobname+"_bcc",atoms=s,concs=c,
                   latpath=latpath,sws=initialsws, xc='PBE')
        swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
        #alloy.lattice_constants_batch_generate(sws=swsrange)
        sws0, B0, e0 = alloy.lattice_constants_analyze(sws=swsrange,prn=False)
        sc_res.append([e0,B0,sws0])
        alloy.bulk(lat='bcc',
                  jobname=finalname+"_bcc",
                   latpath=latpath,
                   sws=sws0,
                   atoms = s,
                   concs = c,
                   splts = split,
                   afm = afm
                   amix=0.02,
                   efmix=0.9,
                   expan='M',
                   sofc='Y',
                   xc='PBE',
                   nky=21)
        alloy.write_inputs()

        # FCC second
        alloy = pyemto.System(folder=apath)
        initialsws = initial_sws[1] # We need some clever way to get this

        alloy.bulk(lat='fcc', jobname=jobname+"_fcc",atoms=s,concs=c,
                   latpath=latpath,sws=initialsws, xc='PBE')
        swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
        
        sws0, B0, e0 = alloy.lattice_constants_analyze(sws=swsrange,prn=False)
        sc_res.append([e0,B0,sws0])
        alloy.bulk(lat='fcc',
                   jobname=finalname+"_fcc",
                   latpath=latpath,
                   sws=sws0,
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
        alloy.write_inputs()

        # HCP last
        alloy = pyemto.System(folder=apath)
        initialsws = initial_sws[2] # We need some clever way to get this
        alloy.bulk(lat='hcp',jobname=jobname,latpath=latpath,
                   sws=initialsws, atoms = s,concs = c, xc='PBE')
        swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
        #alloy.lattice_constants_batch_generate(sws=swsrange)        
        sws0, c_over_a0, B0, e0, R0, cs0 = alloy.lattice_constants_analyze(sws=swsrange,prn=False)
        alloy.sws = sws0 
        ca = round(c_over_a0,3)
        sc_res.append([e0,B0,sws0,c_over_a0])
        # Check is bmdl, kstr and kstr exsist with correct c over a
        hcpname ="hcp_"+str(ca) # Structure name
        strucpath = "../"
        # Check if input files are in place
        if os.path.exists(os.path.join(strucpath,hcpname+".bmdl")):
            pass
        else:
            print("Making structures")
            # make input files
            
            alloy.lattice.set_values(jobname=hcpname,latpath="",
                                         lat='hcp',kappaw=[0.0,-20.0],msgl=0,ca=ca,
                                         dmax=2.2)
            alloy.lattice.bmdl.write_input_file(folder=strucpath)
            alloy.lattice.kstr.write_input_file(folder=strucpath)
            alloy.lattice.shape.write_input_file(folder=strucpath)
            alloy.lattice.batch.write_input_file(folder=strucpath)
        # Make kfcd and kgrn input files
        alloy.bulk(lat='hcp',
                   jobname=finalname+"_hcp",
                   latpath=latpath,
                   latname=hcpname,
                   sws=sws0,
                   ca= ca,
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
        alloy.write_inputs()
    results.append([[s,c],sc_res])

print("Results obtained:")
for r in results:
    # Generate system name
    sname = ""
    for i in range(len(r[0][0])):
        sname=sname+r[0][0][i]+str(r[0][1][i])
    output = "System: "+sname+"\n"
    output = output + "  Magn: " +magn+"\n"
    bcc = r[1][0]
    output = output+"# Strc.     E         sws      B        (c/a)\n"
    output = output+"   bcc: %f %f %f\n" %(bcc[0],bcc[1],bcc[2])
    fcc = r[1][1]
    output = output + "   fcc: %f %f %f\n" %(fcc[0],fcc[1],fcc[2])
    hcp = r[1][2]
    output = output +"   hpc: %f %f %f %f\n" %(hcp[0],hcp[1],hcp[2],hcp[3])
    print(output)
