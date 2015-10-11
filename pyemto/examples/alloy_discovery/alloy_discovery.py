#!/usr/bin/python

import numpy as np
import os

# Help python to find the pyemto folder
#import sys
#sys.path.insert(0, "/home/henrik/local_emto_stuff/pyemto")

import pyemto

##################################################################
##################################################################
##################################################################
# "Input settings" start here.
##################################################################
##################################################################
##################################################################

EMTOdir = "/home/henrik/local_emto_stuff/EMTO6.0"
latpath = "/home/henrik/local_emto_stuff/structures" # Path to bmdl, kstr and shape directories
base_dir = os.getcwd()

# CoCrFeMnNi

# each system need to have same number of alloy elements
systems = [['Co','Cr','Fe','Mn','Ni']]

# Compute concentrations
concs_per_cent = np.array([20.0,15.0,10.0,5.0,0.0])
concentrations = []

for i in range(len(systems[0])):
    for j in range(len(concs_per_cent)):
        if i > 0 and j == 0:
            pass
        else:
            conc = concs_per_cent[j] / 100
            conc_rest = (1.0 - conc) / 4
            if i == 0:
                concentrations.append([conc,conc_rest,conc_rest,conc_rest,conc_rest])
            if i == 1:
                concentrations.append([conc_rest,conc,conc_rest,conc_rest,conc_rest])
            if i == 2:
                concentrations.append([conc_rest,conc_rest,conc,conc_rest,conc_rest])
            if i == 3:
                concentrations.append([conc_rest,conc_rest,conc_rest,conc,conc_rest])
            if i == 4:
                concentrations.append([conc_rest,conc_rest,conc_rest,conc_rest,conc])

# Possible NM (Non-magnetic), FM (ferromagnetic) and 
# DLM (Disordered local moments).
#magn = "NM"
#magn = "FM"
magn = "DLM" 

initial_sws = 2.61

mode = 'create_inputs'
#mode = 'compute_eq_energy'
#mode = 'analyze_results'

##################################################################
##################################################################
##################################################################
# Generic code starts here.
##################################################################
##################################################################
##################################################################

# Check that initialsws is correct format
if type(initial_sws) is float:
    initial_sws = [initial_sws for x in range(3)]
elif type(initial_sws) is list:
    pass
else:
    print("ERROR: Initial_sws should be float or list of 3 floats")
    exit()
if not len(initial_sws) == 3:
    print("ERROR: intial_sws should be a float or list of 3 floats!")
    exit()
            
            
# Sanity checks
for s in systems:
    if not len(s) == len(systems[0]):
        print("Each system need to have same number of alloy elements!")
        exit()
for c in concentrations:
    if not len(c) == len(systems[0]):
        print("Each given concentrations must have same number number as elements in system!")
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
    #for s in systems:
    #    splt = []
    #    for atom in s:
    #        if atom == "Fe":
    #            splt.append(2.0)
    #        else:
    #            splt.append(0.5)
    #    splits.append(splt)
    #
    # Create special DLM lists for atoms and concetration
    #
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
            newc.append(conc/2.0)
            newc.append(conc/2.0)
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# We are ready to make inputs
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if mode == 'create_inputs':
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
        #nlist = s
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
                clist = [2*c[i] for i in range(0,len(c),2)]
            else:
                clist = c
            #clist = c

            for conc in clist:
                count += 1
                cname = cname +str(int(conc*1000)).zfill(4)
                
                if not count == len(clist):
                    cname = cname+"-"

            apath = os.path.join(sname,cname)
            apath = os.path.join(base_dir,apath)
            if not os.path.lexists(apath):
                os.makedirs(apath)
            # Make subdirectory for magnetic state
            apath = os.path.join(apath,magn)
            if not os.path.lexists(apath):
                os.makedirs(apath)

            # Copy QNA parameter file to each folder
            os.system("cp ./QNA_database.dat {0}".format(apath))
            
            # Construct base jobname
            jobname = ""
            
            for i in range(len(nlist)):
                if jobname == "":
                    pass
                else:
                    jobname = jobname + "_"
                jobname = jobname + nlist[i].lower() + "%4.2f" % (clist[i])

            # BCC first
            alloy = pyemto.System(folder=apath,EMTOdir=EMTOdir)
            initialsws = initial_sws[0]
            alloy.bulk(lat='bcc',
                       jobname=jobname+"_bcc",
                       latpath=latpath,
                       msgl=1,
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
                       runtime="12:00:00",
                       stmp='A',
                       niter=300)#,
                       #runKGRN=False)
            swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            alloy.lattice_constants_batch_generate(sws=swsrange)

            # FCC second
            alloy = pyemto.System(folder=apath,EMTOdir=EMTOdir)
            initialsws = initial_sws[1]
            alloy.bulk(lat='fcc',
                       jobname=jobname+"_fcc",
                       latpath=latpath,
                       msgl=1,
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
                       runtime="12:00:00",
                       stmp='A',
                       niter=300)#,
                       #runKGRN=False)
            swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            alloy.lattice_constants_batch_generate(sws=swsrange)

            """
            # HCP last
            alloy = pyemto.System(folder=apath,EMTOdir=EMTOdir)
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
            """

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# We are ready to compute the equilibrium energy
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

elif mode == 'compute_eq_energy':
    results = []
    for si in range(len(systems)):
        s = systems[si]
        split = splits[si]
        # Create main directory 
        sname = ""
        if magn == "DLM":
            nlist = [s[i] for i in range(0,len(s),2)]            
        else:
            nlist = s
        #nlist = s
        for atom in nlist:
            sname = sname + atom
        
        for c in concentrations:
            sc_res = []
            # Make subdirectory for concentration
            cname = ""
            count = 0
            if magn == "DLM":
                clist = [2*c[i] for i in range(0,len(c),2)]
            else:
                clist = c
            #clist = c

            for conc in clist:
                count += 1
                cname = cname +str(int(conc*1000)).zfill(4)
                
                if not count == len(clist):
                    cname = cname+"-"

            apath = os.path.join(sname,cname)
            apath = os.path.join(base_dir,apath)
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
            initialsws = initial_sws[0]
            alloy = pyemto.System(folder=apath,EMTOdir=EMTOdir)
            alloy.bulk(lat='bcc', jobname=jobname+"_bcc",atoms=s,concs=c,
                       latpath=latpath,sws=initialsws, xc='PBE')

            swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            sws0, B0, e0 = alloy.lattice_constants_analyze(sws=swsrange,prn=False)

            alloy.bulk(lat='bcc',
                       jobname=finalname+"_bcc",
                       latpath=latpath,
                       msgl=1,
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
                       runtime="12:00:00",
                       stmp='A',
                       niter=300)#,
                       #runKGRN=False)
            alloy.write_inputs()

            # FCC second
            initialsws = initial_sws[1]
            alloy = pyemto.System(folder=apath,EMTOdir=EMTOdir)
            alloy.bulk(lat='fcc', jobname=jobname+"_fcc",atoms=s,concs=c,
                       latpath=latpath,sws=initialsws, xc='PBE')

            swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            sws0, B0, e0 = alloy.lattice_constants_analyze(sws=swsrange,prn=False)

            alloy.bulk(lat='fcc',
                       jobname=finalname+"_fcc",
                       latpath=latpath,
                       msgl=1,
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
                       runtime="12:00:00",
                       stmp='A',
                       niter=300)#,
                       #runKGRN=False)
            alloy.write_inputs()

            """
            # HCP last
            initialsws = initial_sws[2]
            alloy = pyemto.System(folder=apath)
            initialsws = 3.0 # We need some clever way to get this
            alloy.bulk(lat='hcp',
                       #jobname=jobname+"_hcp",
                       jobname=jobname, # hcp add automatically hcp string to jobname
                       latpath=latpath, sws=initialsws, atoms = s,
                       concs = c, xc='PBE')
     
            swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            #alloy.lattice_constants_batch_generate(sws=swsrange)        
            sws0, c_over_a0, B0, e0, R0, cs0 = alloy.lattice_constants_analyze(sws=swsrange,prn=False)
            alloy.sws = sws0 
            ca = round(c_over_a0,3)
            hcpname ="hcp_"+str(ca) # Structure name
            alloy.bulk(lat='hcp', jobname=finalname+"_hcp",latpath=latpath, latname=hcpname,
                       sws=sws0, ca= ca, atoms = s, concs = c)
            alloy.write_inputs()
            
            # get energy of final
            e_dft = alloy.get_energy()
            sc_res.append([e_dft,sws0,B0,e0,ca])
            """


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# We are ready to analyze the results
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

elif mode == 'analyze_results':
    results = []
    for si in range(len(systems)):
        s = systems[si]
        split = splits[si]
        # Create main directory 
        sname = ""
        if magn == "DLM":
            nlist = [s[i] for i in range(0,len(s),2)]            
        else:
            nlist = s
        #nlist = s
        for atom in nlist:
            sname = sname + atom
        
        for c in concentrations:
            sc_res = []
            # Make subdirectory for concentration
            cname = ""
            count = 0
            if magn == "DLM":
                clist = [2*c[i] for i in range(0,len(c),2)]
            else:
                clist = c
            #clist = c

            for conc in clist:
                count += 1
                cname = cname +str(int(conc*1000)).zfill(4)
                
                if not count == len(clist):
                    cname = cname+"-"

            apath = os.path.join(sname,cname)
            apath = os.path.join(base_dir,apath)
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
            initialsws = initial_sws[0]
            alloy = pyemto.System(folder=apath,EMTOdir=EMTOdir)
            alloy.bulk(lat='bcc', jobname=jobname+"_bcc",atoms=s,concs=c,
                       latpath=latpath,sws=initialsws, xc='PBE')

            swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            sws0, B0, e0 = alloy.lattice_constants_analyze(sws=swsrange,prn=False)

            alloy.bulk(lat='bcc',
                       jobname=finalname+"_bcc",
                       latpath=latpath,
                       sws=sws0,
                       atoms = s,
                       concs = c)

            # get energy of final
            e_dft = alloy.get_energy()
            # Sanity check
            if e_dft == None:
                e_dft = 0.0
            sc_res.append([e_dft,sws0,B0,e0])
     
            # FCC second
            initialsws = initial_sws[1]
            alloy = pyemto.System(folder=apath,EMTOdir=EMTOdir)
            alloy.bulk(lat='fcc', jobname=jobname+"_fcc",atoms=s,concs=c,
                       latpath=latpath,sws=initialsws, xc='PBE')

            swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            sws0, B0, e0 = alloy.lattice_constants_analyze(sws=swsrange,prn=False)

            alloy.bulk(lat='fcc',
                       jobname=finalname+"_fcc",
                       latpath=latpath,
                       sws=sws0,
                       atoms = s,
                       concs = c)

            # get energy of final
            e_dft = alloy.get_energy()
            # Sanity check
            if e_dft == None:
                e_dft = 0.0
            sc_res.append([e_dft,sws0,B0,e0])
     
            """
            # HCP last
            initialsws = initial_sws[2]
            alloy = pyemto.System(folder=apath,EMTOdir=EMTOdir)
            initialsws = 3.0 # We need some clever way to get this
            alloy.bulk(lat='hcp',
                       #jobname=jobname+"_hcp",
                       jobname=jobname, # hcp add automatically hcp string to jobname
                       latpath=latpath, sws=initialsws, atoms = s,
                       concs = c, xc='PBE')
     
            swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            #alloy.lattice_constants_batch_generate(sws=swsrange)        
            sws0, c_over_a0, B0, e0, R0, cs0 = alloy.lattice_constants_analyze(sws=swsrange,prn=False)
            alloy.sws = sws0 
            ca = round(c_over_a0,3)
            hcpname ="hcp_"+str(ca) # Structure name
            alloy.bulk(lat='hcp', jobname=finalname+"_hcp",latpath=latpath, latname=hcpname,
                       sws=sws0, ca= ca, atoms = s, concs = c)
            alloy.write_inputs()
            
            # get energy of final
            e_dft = alloy.get_energy()
            sc_res.append([e_dft,sws0,B0,e0,ca])
            """
            results.append([[s,c],sc_res])


    from pyemto.utilities.utils import wsrad_to_latparam
    
    for r in results:
        # Generate system name
        sname = ""
        for i in range(len(r[0][0])):
            sname=sname+r[0][0][i]+'{0:07.4f}_'.format(100*r[0][1][i])
        # Delete the last unnecessary underscore
        sname = sname[:-1]
        output = "System: "+sname+"\n"
        output = output + "  Magn: " +magn+"\n"
        bcc = r[1][0]
        bcc_lc = wsrad_to_latparam(bcc[1],'bcc')
        output = output+"# Strc.   dft E      lc       sws       B         fit E     (c/a)\n"
        output = output+"   bcc: %f %f %f %f %f %f\n" %(bcc[0],bcc_lc,bcc[1],bcc[2],bcc[3],1.0)
        fcc = r[1][1]
        fcc_lc = wsrad_to_latparam(fcc[1],'fcc')
        output = output + "   fcc: %f %f %f %f %f %f\n" %(fcc[0],fcc_lc,fcc[1],fcc[2],fcc[3],1.0)
        """
        hcp = r[1][2]
        hcp_lc = wsrad_to_latparam(hcp[1],'hcp',ca=hcp[4])
        output = output +"   hpc: %f %f %f %f %f %f\n" %(hcp[0],hcp_lc,hcp[1],hcp[2],hcp[3],hcp[4])
        """
        if magn == "DLM" or magn == "FM":
            # Print magnetic states of system if available
            pass
        
        print(output)
    