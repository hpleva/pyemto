import pyemto
import pyemto.utilities as utils
import numpy as np
import os

latpath = "../../../" # Path do bmdl, kstr and shape directories


# each system need to have same number of alloy elements
#systems = [['Fe','Al'],['Fe','Cr']]
systems = [['Fe'],['Al']]
systems = [['Al']]
#concentrations = [[0.5,0.5]]
concentrations = [[1.0]]

magn = "NM" # Possible values DLM, FM and NM


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
    for s in systems:
        splt = []
        for atom in s:
            if atom == "Fe":
                splt.append(2.0)
            else:
                splt.append(0.5)
        splits.append(splt)
elif magn == "DLM":
    # First duplicate each atoms and concetration
    newsystems = []
    newconcs = []
    for i in range(len(systems)):
        news = []
        newc = []
        splt = []
        for j in range(len(systems[i])):
            print(i,j)
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
else:
    for s in systems:
        splt = []
        for atom in s:
            splt.append(0.0)
        splits.append(splt)


results = []        
#We are ready to make inputs
for si in range(len(systems)):
    s = systems[si]
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
        initialsws = 3.0 # We need some clever way to get this
        alloy.bulk(lat='bcc', jobname=jobname+"_bcc",atoms=s,concs=c,
                   latpath=latpath,sws=initialsws, xc='PBE')
        swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
        #alloy.lattice_constants_batch_generate(sws=swsrange)
        sws0, B0, e0 = alloy.lattice_constants_analyze(sws=swsrange,prn=False)

        alloy.bulk(lat='bcc',jobname=finalname+"_bcc",latpath=latpath,
                   sws=sws0,atoms = s,concs = c)

        # get energy of final
        e_dft = alloy.get_energy()
        sc_res.append([e_dft,sws0,B0,e0])


        # FCC second
        alloy = pyemto.System(folder=apath)
        initialsws = 3.0 # We need some clever way to get this

        alloy.bulk(lat='fcc', jobname=jobname+"_fcc",atoms=s,concs=c,
                   latpath=latpath,sws=initialsws, xc='PBE')
        swsrange = np.linspace(initialsws-0.1,initialsws+0.1,7) # A list of 7 different volumes
        
        sws0, B0, e0 = alloy.lattice_constants_analyze(sws=swsrange,prn=False)
        alloy.bulk(lat='fcc', jobname=finalname+"_fcc", latpath=latpath,sws=sws0,
                   atoms = s, concs = c)

        # get energy of final
        e_dft = alloy.get_energy()
        sc_res.append([e_dft,sws0,B0,e0])

        
        # HCP last
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
    results.append([[s,c],sc_res])
    
for r in results:
    # Generate system name
    sname = ""
    for i in range(len(r[0][0])):
        sname=sname+r[0][0][i]+str(r[0][1][i])
    output = "System: "+sname+"\n"
    output = output + "  Magn: " +magn+"\n"
    bcc = r[1][0]
    bcc_lc = utils.wsrad_to_latparam(bcc[1],'bcc')
    output = output+"# Strc.   dft E      lc       sws       B         fit E     (c/a)\n"
    output = output+"   bcc: %f %f %f %f %f\n" %(bcc[0],bcc_lc,bcc[1],bcc[2],bcc[3])
    fcc = r[1][1]
    fcc_lc = utils.wsrad_to_latparam(fcc[1],'fcc')
    output = output + "   fcc: %f %f %f %f %f\n" %(fcc[0],fcc_lc,fcc[1],fcc[2],fcc[3])
    hcp = r[1][2]
    hcp_lc = utils.wsrad_to_latparam(hcp[1],'hcp',ca=hcp[4])
    output = output +"   hpc: %f %f %f %f %f %f\n" %(hcp[0],hcp_lc,hcp[1],hcp[2],hcp[3],hcp[4])

    if magn == "DLM" or magn == "FM":
        # Print magnetic states of system if available
        pass
    
    
    print(output)

