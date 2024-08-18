import flopy
import numpy as np
import os
import post_processing as proc

def create_transport_package(name, gwt, type, transport_kwargs, strt_cim=0):

    if type=='dispersive':
        flopy.mf6.ModflowGwtdsp(gwt,
                                xt3d_off=True,
                                diffc=8.64e-5,
                                alh=transport_kwargs['alpha'],
                                ath1=transport_kwargs['alpha'] * 0.1,
                                alv=transport_kwargs['alpha'] * 0.01)

    elif type=='dual porosity':
        flopy.mf6.ModflowGwtdsp(gwt,
                                xt3d_off=True,
                                diffc=8.64e-5,
                                alh=transport_kwargs['alpha'],
                                ath1=transport_kwargs['alpha'] * 0.1,
                                alv=transport_kwargs['alpha'] * 0.01)

        flopy.mf6.ModflowGwtist(gwt,
                                cim_filerecord=f"{name}.cim",
                                porosity=transport_kwargs['porosity_im'],
                                thetaim=transport_kwargs['porosity_im'],
                                volfrac=transport_kwargs['volfrac'],
                                zetaim=transport_kwargs['zetaim'],
                                cim=35
        )

def create_k_data(H, D, ncol, nlay, ka, kb, anis, ka_anis=1):
    hk = np.concatenate((kb * np.ones((int(nlay*D/(D+H)), 1, ncol)),
                         ka * np.ones((int(nlay*H/(D+H)), 1, ncol))), axis=0)
    hka = [anis for i in range(int(nlay*D/(D+H)))] + [ka_anis for i in range(int(nlay*H/(D+H)))]
    return hk, hka

def build_2D_steady_model(dir, name, type, base_params, transport_kwargs):
    D = 20
    H = 20
    L = 4400
    delc = 12.5
    delr = 12.5
    delv = 2.5
    nlay = int((D + H) / delv)
    ncol = int(L / delc)
    nrow = 1
    perlen = 1e8
    nper = 1
    nstp = 1e5
    tsmult = 1
    x_onshore = 400
    h_onshore = 0.75
    anis = 0.01
    kb = 0.1

    mf6exe = "/home/superuser/mf6"
    # mf6exe = "C:/Users/ccl124/code/code_cdrive/bin/mf6.exe"
    length_units = "m"
    time_units = "days"

    nouter, ninner = 500, 500
    hclose, rclose, relax = 1e-5, 1e-5, 0.97

    modelname = name

    ws = str(dir)
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
    botm = np.linspace(-delv, -int(nlay * delv), nlay)
    idomain = np.ones((nlay, nrow, ncol), dtype=np.int32)
    idomain[int(D / delv), :, 0] = 1
    idomain[:, :, -1] = 1
    idomain[0, 0, int(x_onshore / delc):] = 1

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
    hk, hka = create_k_data(H, D, ncol, nlay, base_params['ka'], kb, anis)
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

    flopy.mf6.ModflowGwfic(gwf, strt=np.ones_like(idomain) * 0.0)

    pd = [(0, 0.7, 0.0, "trans", "concentration")]

    flopy.mf6.ModflowGwfbuy(gwf, packagedata=pd)

    ghb_data = {}
    chd_data = {}
    ghb_temp = []
    chd_temp = []

    # add onshore boundary
    for lay in range(int(D / delv), int((D + H) / delv)):
        for row in range(nrow):
            chd_temp.append([(lay, row, 0), h_onshore, 0.0])
    # add vertical
    for lay in range(int((D + H) / delv)):
        for row in range(nrow):
            cond = hk[lay, row, ncol - 1] * delv * delc / (delr)
            ghb_temp.append([(lay, row, ncol - 1), 0, cond, 35.0])
    # add horizontal
    for col in range(ncol - int(x_onshore / delc), ncol):
        for row in range(nrow):
            cond = hk[0, row, col] * hka[0] * delr * delc / (delv)
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

    strt = np.ones_like(idomain) * 35.0
    strt[:, :, :int(x_onshore / L * ncol)] = 0.0
    flopy.mf6.ModflowGwtic(gwt, strt=strt)

    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    create_transport_package(name, gwt, type, transport_kwargs)

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

def build_2D_transient_model(name, conc, head, type, base_params, transport_kwargs):
    pass

def get_steady_metrics_and_arrays(dir, name):

    try:
        sim = flopy.mf6.MFSimulation.load(
            sim_ws=str(dir),
            exe_name="/home/superuser/mf6",
            verbosity_level=0,
        )
        gwf = sim.get_model(name)
        gwt = sim.get_model("trans")
        bud = gwf.output.budget()
        times = bud.get_times()
        time = times[-1]
        conc = gwt.output.concentration().get_data(totim=time)
        toe_mean = proc.find_average_toe_position(conc, 12.5, 400)
        width_mean = proc.find_average_width(conc, 12.5, int(352), 400)
        head = gwt.output.head().get_data(totim=time)
        return toe_mean, width_mean, conc, head

    except:
        return np.inf, np.inf, np.inf, np.inf