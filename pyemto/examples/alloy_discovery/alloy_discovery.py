#!/usr/bin/python

import numpy as np
import os

# Help python to find the pyemto folder
#import sys
#sys.path.insert(0, "/home/hpleva/pyemto")
import pyemto


##################################################################
##################################################################
##################################################################
# "Input settings" start here.
##################################################################
##################################################################
##################################################################

EMTOdir = "/home/hpleva/EMTO/openmp-stable-cmake"
latpath = "/home/hpleva/structures" # Path to bmdl, kstr and shape directories
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
            conc = concs_per_cent[j] #/ 100
            conc_rest = (100.0 - conc) / 4
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

# Pure elements
concentrations.append([100.0, 0, 0, 0, 0])
concentrations.append([0, 100.0, 0, 0, 0])
concentrations.append([0, 0, 100.0, 0, 0])
concentrations.append([0, 0, 0, 100.0, 0])
concentrations.append([0, 0, 0, 0, 100.0])

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
                       msgl=0,
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
                       niter=300,
                       tole=1.0e-8,
                       nky=31,
                       alpcpa=0.796)#,
                       #runKGRN=False)
            #swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            swsrange = np.linspace(2.40,2.80,20)
            alloy.lattice_constants_batch_generate(sws=swsrange)

            # FCC second
            alloy = pyemto.System(folder=apath,EMTOdir=EMTOdir)
            initialsws = initial_sws[1]
            alloy.bulk(lat='fcc',
                       jobname=jobname+"_fcc",
                       latpath=latpath,
                       msgl=0,
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
                       niter=300,
                       nky=25,
                       tole=1.0e-8,
                       alpcpa=0.796)#,
                       #runKGRN=False)
            #swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            # FM
            #swsrange = np.linspace(2.56,2.66,15)
            #swsrange = np.linspace(2.56,2.76,29)
            # DLM
            swsrange = np.linspace(2.40,2.80,20)
            alloy.lattice_constants_batch_generate(sws=swsrange)

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
                       nkz=13,
                       stmp='A',
                       niter=300,
                       runtime='24:00:00',
                       tole=1.0e-8,
                       alpcpa=0.796)#,
            #swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            swsrange = np.linspace(2.40,2.80,20)
            alloy.lattice_constants_batch_generate(sws=swsrange)

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
                       msgl=0,
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
                       niter=300,
                       alpcpa=0.796)#,
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
                       msgl=0,
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
                       niter=300,
                       alpcpa=0.796)#,
                       #runKGRN=False)
            alloy.write_inputs()

            # HCP last
            initialsws = initial_sws[2]
            alloy = pyemto.System(folder=apath)
            alloy.bulk(lat='hcp',
                       #jobname=jobname+"_hcp",
                       jobname=jobname, # hcp add automatically hcp string to jobname
                       latpath=latpath,
                       sws=initialsws,
                       atoms = s,
                       concs = c, 
                       xc='PBE')
     
            swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            #alloy.lattice_constants_batch_generate(sws=swsrange)        
            sws0, c_over_a0, B0, e0, R0, cs0 = alloy.lattice_constants_analyze(sws=swsrange,prn=False)
            alloy.sws = sws0 
            ca = round(c_over_a0,3)
            hcpname ="hcp_"+str(ca) # Structure name
            structpath = "../"
            # Check if input files are in place
            if os.path.exists(os.path.join(structpath,hcpname+".bmdl")):
                pass
            else:
                print("Making structures")
                # make input files
                alloy.lattice.set_values(jobname=hcpname,latpath="",
                                         lat='hcp',kappaw=[0.0,-20.0],msgl=0,ca=ca,
                                         dmax=2.2)
                alloy.lattice.bmdl.write_input_file(folder=structpath)
                alloy.lattice.kstr.write_input_file(folder=structpath)
                alloy.lattice.shape.write_input_file(folder=structpath)
                alloy.lattice.batch.write_input_file(folder=structpath)
                

            # Make kfcd and kgrn input files
            print("hcp", afm)
            alloy.bulk(lat='hcp',
                       #jobname=jobname+"_hcp",
                       jobname=finalname+"_hcp",
                       latpath=latpath,
                       latname=hcpname,
                       sws=sws0,
                       ca = ca,
                       msgl = 0,
                       atoms = s,
                       concs = c,
                       splts = split,
                       afm = afm,
                       amix=0.02,
                       efmix=0.9,
                       expan='M',
                       sofc='Y',
                       xc='PBE',
                       nky=11,
                       nkz=7,
                       alpcpa=0.796)

            alloy.write_inputs()
            
           

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

            """
            # BCC first
            initialsws = initial_sws[0]
            alloy = pyemto.System(folder=apath,EMTOdir=EMTOdir)
            alloy.bulk(lat='bcc', jobname=jobname+"_bcc",atoms=s,concs=c,
                       latpath=latpath,sws=initialsws, xc='PBE')

            swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            sws0, B0, e0, grun0, R_squared = alloy.lattice_constants_analyze(sws=swsrange,prn=False,return_error=True)

            alloy.bulk(lat='bcc',
                       jobname=finalname+"_bcc",
                       latpath=latpath,
                       sws=sws0,
                       atoms = s,
                       concs = c)

            #from pyemto.utilities.utils import run_bash
            #tmp = run_bash("mv "+apath+"/kfcd/"+finalname+"_bcc*.prn"+" "+apath+"/kfcd/"+finalname+"_bcc_"+"{0:8.6f}.prn".format(sws0))
            #tmp = run_bash("mv "+apath+"/kgrn/"+finalname+"_bcc*.prn"+" "+apath+"/kgrn/"+finalname+"_bcc_"+"{0:8.6f}.prn".format(sws0))
            # get energy of final
            e_dft = alloy.get_energy()
            # get magnetic moments
            magn_moms = alloy.get_moments()
            # get total DOS at Fermi level
            dos_total = alloy.get_fdos()            
            # Sanity checks
            if e_dft == None:
                e_dft = 0.0
            if magn_moms == None:
                magn_moms = [0 for i in range(len(concs))]
            sc_res.append([e_dft,sws0,B0,grun0,dos_total,e0,R_squared,magn_moms])
     
            # FCC second
            initialsws = initial_sws[1]
            alloy = pyemto.System(folder=apath,EMTOdir=EMTOdir)
            alloy.bulk(lat='fcc', jobname=jobname+"_fcc",atoms=s,concs=c,
                       latpath=latpath,sws=initialsws, xc='PBE')

            swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            sws0, B0, e0, grun0, R_squared = alloy.lattice_constants_analyze(sws=swsrange,prn=False,return_error=True)

            alloy.bulk(lat='fcc',
                       jobname=finalname+"_fcc",
                       latpath=latpath,
                       sws=sws0,
                       atoms = s,
                       concs = c)

            #tmp = run_bash("mv "+apath+"/kfcd/"+finalname+"_fcc*.prn"+" "+apath+"/kfcd/"+finalname+"_fcc_"+"{0:8.6f}.prn".format(sws0))
            #tmp = run_bash("mv "+apath+"/kgrn/"+finalname+"_fcc*.prn"+" "+apath+"/kgrn/"+finalname+"_fcc_"+"{0:8.6f}.prn".format(sws0))
            # get energy of final
            e_dft = alloy.get_energy()
            # get magnetic moments
            magn_moms = alloy.get_moments()
            # get total DOS at Fermi level
            dos_total = alloy.get_fdos()
            # Sanity check
            if e_dft == None:
                e_dft = 0.0
            if magn_moms == None:
                magn_moms = [0 for i in range(len(concs))]
            sc_res.append([e_dft,sws0,B0,grun0,dos_total,e0,R_squared,magn_moms])
            """
            
            # HCP last
            initialsws = initial_sws[2]
            alloy = pyemto.System(folder=apath,EMTOdir=EMTOdir)
            alloy.bulk(lat='hcp',
                       #jobname=jobname+"_hcp",
                       jobname=jobname, # hcp add automatically hcp string to jobname
                       latpath=latpath, sws=initialsws, atoms = s,
                       concs = c, xc='PBE')
     
            swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
            #alloy.lattice_constants_batch_generate(sws=swsrange)        
            sws0, c_over_a0, B0, e0, R0, cs0, grun0, R_squared = alloy.lattice_constants_analyze(sws=swsrange,prn=True,return_error=True)
            alloy.sws = sws0 
            ca = round(c_over_a0,3)
            hcpname ="hcp_"+str(ca) # Structure name
            alloy.bulk(lat='hcp', jobname=finalname+"_hcp",latpath=latpath, latname=hcpname,
                       sws=sws0, ca= ca, atoms = s, concs = c)

            #tmp = run_bash("mv "+apath+"/kfcd/"+finalname+"_hcp*.prn"+" "+apath+"/kfcd/"+finalname+"_hcp_"+"{0:8.6f}.prn".format(sws0))
            #tmp = run_bash("mv "+apath+"/kgrn/"+finalname+"_hcp*.prn"+" "+apath+"/kgrn/"+finalname+"_hcp_"+"{0:8.6f}.prn".format(sws0))
            # get energy of final
            e_dft = alloy.get_energy()
            # get magnetic moments
            magn_moms = alloy.get_moments()
            # get total DOS at Fermi level
            dos_total = alloy.get_fdos()
            # Sanity check
            if e_dft == None:
                e_dft = 0.0
            if magn_moms == None:
                magn_moms = [0 for i in range(len(concs))]
            sc_res.append([e_dft,sws0,B0,grun0,dos_total,e0,R_squared,magn_moms,ca])

            results.append([[s,c],sc_res])


    from pyemto.utilities.utils import wsrad_to_latparam
    
    output_all = ""

    for r in results:
        # Generate system name and concentrations line
        sname = ""
        conc_line = ""
        for i in range(len(r[0][0])):
            sname = sname + r[0][0][i] + '{0:07.4f}_'.format(100*r[0][1][i])
            conc_line = conc_line + '{0:8.6f} '.format(r[0][1][i])
        # Delete the last unnecessary underscore
        sname = sname[:-1]
        output = "  System: "+sname+"\n"
        output = output + "    Magn: " + magn + "\n"
        output = output + "    Conc: " + conc_line + "\n"
        
        bcc = r[1][0]
        bcc_lc = wsrad_to_latparam(bcc[1],'bcc')
        output = output + "#  Struc.   dft E(Ry)    lc(AA)    sws(bohr) B(GPa)     grun       DOSEf(1/eV)    fit E(Ry)   fit err      (c/a)\n"
        output = output + "     bcc: {0:13.6f} {1:9.6f} {2:9.6f} {3:11.6f} {4:8.6f} {5:10.6f} {6:13.6f} {7:12.10f} {8:8.6f}\n".format(bcc[0],bcc_lc,bcc[1],bcc[2],bcc[3],bcc[4],bcc[5],bcc[6],1.0)
        # Generate the output line for magnetic moments
        output = output + "bcc_moms:"
        for i in range(len(bcc[7])):
            output = output + " {0:9.6f}".format(bcc[7][i])
        output = output + "\n"
        
        fcc = r[1][1]
        fcc_lc = wsrad_to_latparam(fcc[1],'fcc')
        output = output + "     fcc: {0:13.6f} {1:9.6f} {2:9.6f} {3:11.6f} {4:8.6f} {5:10.6f} {6:13.6f} {7:12.10f} {8:8.6f}\n".format(fcc[0],fcc_lc,fcc[1],fcc[2],fcc[3],fcc[4],fcc[5],fcc[6],1.0)
        # Generate the output line for magnetic moments
        output = output + "fcc_moms:"
        for i in range(len(fcc[7])):
            output = output + " {0:9.6f}".format(fcc[7][i])
        output = output + "\n"
        
        hcp = r[1][2]
        hcp_lc = wsrad_to_latparam(hcp[1],'hcp',ca=hcp[8])
        output = output +"     hpc: {0:13.6f} {1:9.6f} {2:9.6f} {3:11.6f} {4:8.6f} {5:10.6f} {6:13.6f} {7:12.10f} {8:8.6f}\n".format(hcp[0],hcp_lc,hcp[1],hcp[2],hcp[3],hcp[4],hcp[5],hcp[6],hcp[8])
        # Generate the output line for magnetic moments
        output = output + "hcp_moms:"
        for i in range(len(hcp[7])):
            output = output + " {0:9.6f}".format(hcp[7][i])
        output = output + "\n"

        output_all += output+"\n"
    output_file = open('results','w')
    output_file.write(output_all)
    output_file.close()
    
