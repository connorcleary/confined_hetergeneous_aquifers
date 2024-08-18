import os
import flopy
import numpy as np
import random_fields as rf

def build_base_model(name, field_name):
    # create model to run to steady state: a semiconfined coastal aquifer with modern OFG
    
    D = 20
    H = 20
    W = 1000
    L = 4400
    delc = 12.5
    delv = 2.5
    delr = 12.5
    nlay = int((D+H)/delv)
    ncol = int(L/delc)
    nrow = int(W/delr)
    perlen = 1e8
    nper = 1
    nstp = 1e5
    tsmult = 1
    x_onshore = 400
    h_onshore = 0.75
    anis = 0.01
    alpha = 0
    kb = 0.1

    mf6exe = "/home/superuser/mf6"
    # mf6exe = "C:/Users/ccl124/code/code_cdrive/bin/mf6.exe"
    length_units = "m"
    time_units = "days"

    nouter, ninner = 500, 500
    hclose, rclose, relax = 1e-5, 1e-5, 0.97

    modelname = name

    ws = f"./model_files/{name}"
    if not os.path.exists(ws):
        os.makedirs(ws)

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name=mf6exe
    )

    time_spd = [(perlen, nstp, tsmult)]

    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=time_spd, time_units=time_units
    )

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwf.name),
    )

    sim.register_ims_package(ims, [gwf.name])

    top = 0
    botm = np.linspace(-delv, -int(nlay*delv), nlay)
    idomain = np.ones((nlay, nrow, ncol), dtype=np.int32)
    idomain[int(D/delv), :, 0] = 1
    idomain[:, :, -1] = 1
    idomain[0, 0, int(x_onshore/delc):] = 1


    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain
    )
    hk, hka, = rf.import_field(field_name, H, delv, D, kb, W, delr, delc=delc)

    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=1,
        k=hk,
        k33overk=True,
        k33=hka,
        rewet_record="REWET WETFCT 1.0 IWETIT 1 IHDWET 0",
        wetdry=1
    )


    flopy.mf6.ModflowGwfic(gwf, strt=np.ones_like(idomain)*0.0)

    pd = [(0, 0.7, 0.0, "trans", "concentration")]

    flopy.mf6.ModflowGwfbuy(gwf, packagedata=pd)

    ghb_data = {}
    chd_data= {}
    ghb_temp = []
    chd_temp = []

    # add onshore boundary
    for lay in range(int(D/delv), int((D+H)/delv)):
        for row in range(nrow):
            chd_temp.append([(lay, row, 0), h_onshore, 0.0])
    # add vertical
    for lay in range(int((D+H)/delv)):
        for row in range(nrow):
            cond = hk[lay, row, ncol-1] * delv * delc / (delr)
            ghb_temp.append([(lay, row, ncol-1), 0, cond, 35.0])
    # add horizontal
    for col in range(ncol-int(x_onshore/delc), ncol):
        for row in range(nrow):
            cond = hk[0, row, col]*hka[0]* delr * delc / (delv)
            ghb_temp.append([(0, row, col), 0, cond, 35.0])

    ghb_data[0] = ghb_temp
    chd_data[0] = chd_temp

    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chd_data,
        pname="CHD-1",
        auxiliary="CONCENTRATION",
    )

    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghb_data,
        pname="GHB-1",
        auxiliary=["CONCENTRATION"]
    )

    ws_results = f"./model_files/{name}"
    if not os.path.exists(ws_results):
        os.makedirs(ws_results)

    head_filerecord = f"{name}.hds"
    budget_filerecord = f"{name}.bud"
    saverecord = {0: [("HEAD", "LAST"), ("BUDGET", "LAST")]}
    saverecord = {0: [("HEAD", "LAST"), ("BUDGET", "LAST")]}

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=saverecord,
    )

    flopy.mf6.ModflowGwfsto(gwf, ss=1e-6, sy=0.3, iconvert=1, transient=True)
    gwt = flopy.mf6.ModflowGwt(sim, modelname="trans")

    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwt.name),
    )

    sim.register_ims_package(imsgwt, [gwt.name])

    flopy.mf6.ModflowGwtdis(
        gwt,
        idomain=idomain,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    flopy.mf6.ModflowGwtmst(gwt, porosity=0.5)

    strt = np.ones_like(idomain)*35.0
    strt[:, :, :int(x_onshore/L*ncol)] = 0.0
    flopy.mf6.ModflowGwtic(gwt, strt=strt)

    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, diffc=8.64e-5, alh=alpha, ath1=alpha*0.1, alv=alpha*0.01)

    sourcerecarray = [
        ("GHB-1", "AUX", "CONCENTRATION"),
        ("CHD-1", "AUX", "CONCENTRATION")
    ]

    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)
    saverecord = {0: [("CONCENTRATION", "LAST")]}

    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{name}.cbc",
        concentration_filerecord=f"{name}.ucn",
        concentrationprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=saverecord,
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=gwt.name
    )

    return sim

def run_model(swt):
    swt.write_simulation(silent=True)
    success, buff = swt.run_simulation(silent=True)
    if not success:
            print(buff)

    return success

def add_pumping_to_model():
    pass

def add_head_to_model(name, field_name, conc, head):

    H = 20
    D = 20
    W = 1000
    L = 4400
    delc = 12.5
    delv = 2.5
    delr = 12.5
    nlay = int((D+H)/delv)
    ncol = int(L/delc)
    nrow = int(W/delr)
    perlen = 365*100
    nper = 1
    nstp = 3650
    tsmult = 1
    x_onshore = 400
    h_onshore = 0
    anis = 0.01
    alpha = 0
    kb = 0.1
    h_modern = -1
    mf6exe = "/home/superuser/mf6"
    # mf6exe = "C:/Users/ccl124/code/code_cdrive/bin/mf6.exe"
    length_units = "m"
    time_units = "days"

    nouter, ninner = 500, 500
    hclose, rclose, relax = 1e-5, 1e-5, 0.97

    modelname = name 

    ws = f"./model_files/{name}"
    if not os.path.exists(ws):
        os.makedirs(ws)

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name=mf6exe
    )

    time_spd = [(perlen, nstp, tsmult)]

    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=time_spd, time_units=time_units
    )

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwf.name),
    )

    sim.register_ims_package(ims, [gwf.name])

    top = 0
    botm = np.linspace(-delv, -int(nlay*delv), nlay)
    idomain = np.ones((nlay, nrow, ncol), dtype=np.int32)
    idomain[int(D/delv), :, 0] = 1
    idomain[:, :, -1] = 1
    idomain[0, 0, int(x_onshore/delc):] = 1

    
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain
    )    
    hk, hka, = rf.import_field(field_name, H, delv, D, kb, W, delr, delc=delc)

    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=1,
        k=hk,
        k33overk=True, 
        k33=hka,
        rewet_record="REWET WETFCT 1.0 IWETIT 1 IHDWET 0",
        wetdry=1
    )


    flopy.mf6.ModflowGwfic(gwf, strt=head)

    pd = [(0, 0.7, 0.0, "trans", "concentration")]

    flopy.mf6.ModflowGwfbuy(gwf, packagedata=pd)

    ghb_data = {}
    chd_data= {}
    ghb_temp = []
    chd_temp = []

    # add onshore boundary 
    for lay in range(int(D/delv), int((D+H)/delv)):
        for row in range(nrow):
            chd_temp.append([(lay, row, 0), h_modern, 0.0])
    # add vertical 
    for lay in range(int((D+H)/delv)): 
        for row in range(nrow):
            cond = hk[lay, row, ncol-1] * delv * delc / (delr)
            ghb_temp.append([(lay, row, ncol-1), 0, cond, 35.0])
    # add horizontal 
    for col in range(ncol-int(x_onshore/delc), ncol):
        for row in range(nrow):
            cond = hk[0, row, col]*hka[0]* delr * delc / (delv)
            ghb_temp.append([(0, row, col), 0, cond, 35.0])

    ghb_data[0] = ghb_temp
    chd_data[0] = chd_temp

    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chd_data,
        pname="CHD-1",
        auxiliary="CONCENTRATION",
    )

    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghb_data,
        pname="GHB-1",
        auxiliary=["CONCENTRATION"]
    )
    
    ws_results = f"./model_files/{name}"
    if not os.path.exists(ws_results):
        os.makedirs(ws_results)

    head_filerecord = f"{name}.hds"
    budget_filerecord = f"{name}.bud"
    saverecord = {0: [("HEAD", "LAST"), ("BUDGET", "LAST")]}
    saverecord = {0: [("HEAD", "LAST"), ("BUDGET", "LAST")]}

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=saverecord,
    )

    flopy.mf6.ModflowGwfsto(gwf, ss=1e-6, sy=0.3, iconvert=1, transient=True)
    gwt = flopy.mf6.ModflowGwt(sim, modelname="trans")
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwt.name),
    )

    sim.register_ims_package(imsgwt, [gwt.name])

    flopy.mf6.ModflowGwtdis(
        gwt,
        idomain=idomain,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    flopy.mf6.ModflowGwtmst(gwt, porosity=0.5)

    flopy.mf6.ModflowGwtic(gwt, strt=conc)

    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, diffc=8.64e-5, alh=alpha, ath1=alpha*0.1, alv=alpha*0.01)

    sourcerecarray = [
        ("GHB-1", "AUX", "CONCENTRATION"),
        ("CHD-1", "AUX", "CONCENTRATION")
    ]

    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)
    # save result every 10 years
    saverecord = {0: [("CONCENTRATION", "FREQUENCY", 365)]}
    	
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{name}.cbc",
        concentration_filerecord=f"{name}.ucn",
        concentrationprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=saverecord,
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=gwt.name
    )

    return sim

def add_recovery_to_model():
    pass

def get_bulk_conductivities(L, H, W, i):

    name = "temp"
    mf6exe = "/home/superuser/mf6"
    # mf6exe = "C:/Users/ccl124/code/code_cdrive/bin/mf6.exe"
    length_units = "m"
    time_units = "days"

    delc = 6.25
    delv = 2.5
    delr = 6.25
    nlay = int(H/delv)
    ncol = int(L/delc)
    nrow = int(W/delr)
    nouter, ninner = 500, 500
    hclose, rclose, relax = 1e-5, 1e-5, 0.97

    modelname = name 

    ws = f"./model_files/{name}"
    if not os.path.exists(ws):
        os.makedirs(ws)

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name=mf6exe
    )

    time_spd = [(1.0, 1, 1.0)]

    flopy.mf6.ModflowTdis(
        sim, nper=1, perioddata=time_spd, time_units=time_units
    )

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True, print_flows=True)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwf.name),
    )

    sim.register_ims_package(ims, [gwf.name])

    top = 0
    botm = np.linspace(-delv, -int(nlay*delv), nlay)
    idomain = np.ones((nlay, nrow, ncol), dtype=np.int32)


    
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain
    )    
    # field_name, H, delv, D, kb, W=3600, delr=12.5, L=4400, delc=12.5
    hk, hka, = rf.import_bulk_field(f"TSim_Out{i+1}", H, delv, 0, 0, W, delr, L, delc=delc)

    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=1,
        k=hk,
        k33overk=True, 
        k33=hka,
        rewet_record="REWET WETFCT 1.0 IWETIT 1 IHDWET 0",
        wetdry=1
    )


    flopy.mf6.ModflowGwfic(gwf, strt=np.ones_like(idomain)*0.0)

    chd_data= {}
    chd_temp = []

    # add onshore boundary 
    for lay in range(int(H/delv)):
        for row in range(nrow):
            chd_temp.append([(lay, row, 0), 1])
    # add vertical 
    for lay in range(int((H)/delv)): 
        for row in range(nrow):
            chd_temp.append([(lay, row, ncol-1), 0])

    chd_data[0] = chd_temp

    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chd_data,
        pname="CHD-1",
    )

    ws_results = f"./model_files/{name}"
    if not os.path.exists(ws_results):
        os.makedirs(ws_results)

    head_filerecord = f"{name}.hds"
    budget_filerecord = f"{name}.bud"
    saverecord = {0: [("HEAD", "LAST"), ("BUDGET", "LAST")]}

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=saverecord,
    )

    sim.write_simulation()
    success, buff = sim.run_simulation(silent=False)

    bud = gwf.output.budget()
    times = bud.get_times()
    time = times[-1]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
                bud.get_data(text="DATA-SPDIS", totim=time)[0],
                gwf,
               )
    
    outflow = np.sum(qx[:, :, -1])*delr*delv
    hk_eff = outflow/W/H/(1/L)

    
    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name=mf6exe
    )

    time_spd = [(1.0, 1, 1.0)]

    flopy.mf6.ModflowTdis(
        sim, nper=1, perioddata=time_spd, time_units=time_units
    )

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True, print_flows=True)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwf.name),
    )

    sim.register_ims_package(ims, [gwf.name])

    top = 0
    botm = np.linspace(-delv, -int(nlay*delv), nlay)
    idomain = np.ones((nlay, nrow, ncol), dtype=np.int32)


    
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain
    )    
    # field_name, H, delv, D, kb, W=3600, delr=12.5, L=4400, delc=12.5
    hk, hka, = rf.import_bulk_field(f"TSim_Out{i+1}", H, delv, 0, 0, W, delr, L, delc=delc)

    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=1,
        k=hk,
        k33overk=True, 
        k33=hka,
        rewet_record="REWET WETFCT 1.0 IWETIT 1 IHDWET 0",
        wetdry=1
    )


    flopy.mf6.ModflowGwfic(gwf, strt=np.ones_like(idomain)*0.0)

    chd_data= {}
    chd_temp = []

    # add above boundary 
    for col in range(int(L/delc)):
        for row in range(nrow):
            chd_temp.append([(0, row, col), 0])
    # add below
    for col in range(int(L/delc)):
        for row in range(nrow):
            chd_temp.append([(nlay-1, row, col), 1])

    chd_data[0] = chd_temp

    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chd_data,
        pname="CHD-1",
    )

    ws_results = f"./model_files/{name}"
    if not os.path.exists(ws_results):
        os.makedirs(ws_results)

    head_filerecord = f"{name}.hds"
    budget_filerecord = f"{name}.bud"
    saverecord = {0: [("HEAD", "LAST"), ("BUDGET", "LAST")]}

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=saverecord,
    )

    sim.write_simulation()
    success, buff = sim.run_simulation(silent=False)

    bud = gwf.output.budget()
    times = bud.get_times()
    time = times[-1]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
                bud.get_data(text="DATA-SPDIS", totim=time)[0],
                gwf,
               )
    
    outflow = np.sum(qz[0, :, :])*delr*delc
    hka_eff = outflow/W/L/(1/H)

    return hk_eff, hka_eff


def get_last_results(name):
    ws = f"model_files/{name}"
    if not os.path.exists(ws):
        os.makedirs(ws)


    # .set_trace()
    try:
        sim = flopy.mf6.MFSimulation.load(
            sim_ws=ws,
            exe_name="/home/superuser/mf6",
            verbosity_level=0,
        )
        gwf = sim.get_model(name)
        gwt = sim.get_model("trans")
        bud = gwf.output.budget()                
        times = bud.get_times()
        time = times[-1]
        conc = gwt.output.concentration().get_data(totim=time)
        head = gwt.output.head().get_data(totim=time)
        qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
                    bud.get_data(text="DATA-SPDIS", totim=time)[0],
                    gwf,
                )            
        return conc, head, qx, qy, qz
    
    except:
        return None, None, None, None, None

def get_all_conc(name):

    ws = f"model_files/{name}"
    sim = flopy.mf6.MFSimulation.load(
        sim_ws=ws,
        exe_name="/home/superuser/mf6",
        verbosity_level=0,
    )
    gwf = sim.get_model(name)
    gwt = sim.get_model("trans")
    conc = []
    for time in gwt.output.concentration().get_times():
        conc.append(gwt.output.concentration().get_data(totim=time))
  
    return conc

def run_modpath_modern(params):

    name, interface = params
    porosity = 0.5
    ws = f"model_files/{name}"
    sim = flopy.mf6.MFSimulation.load(
        sim_ws=ws,
        exe_name="/home/superuser/mf6",
        verbosity_level=0,
    )
    gwf = sim.get_model(name)


    mp = flopy.modpath.Modpath7(
        modelname=f"{name}",
        flowmodel=gwf,
        exe_name="/home/superuser/mp7",
        model_ws=ws,
    )
    mpbas = flopy.modpath.Modpath7Bas(mp, porosity=0.5)

    particle_list = []
    for row, col in enumerate(interface):
        for lay in range(8, 16):
            particle_list.append((lay, row, col))

    mp7_particle_data = flopy.modpath.ParticleData(
        particle_list,
        structured=True,
        particleids=[i for i in range(len(particle_list))],
    )

    pg = flopy.modpath.ParticleGroup(
        particlegroupname="PG",
        particledata=mp7_particle_data,
        filename=f"{name}" + "pg.sloc")

    mpsim = flopy.modpath.Modpath7Sim(
        mp,
        pathlinefilename=f"{name}.path",
        simulationtype="pathline",
        trackingdirection="forward",
        weaksinkoption="pass_through",
        weaksourceoption="pass_through",
        budgetoutputoption="summary",
        particlegroups=[pg],
        stoptimeoption='total'
    )

    mp.write_input()
    mp.run_model(silent=False, report=True)

    print("all good")
