import matplotlib.pyplot as plt
import numpy as np
import modelling
import post_processing as proc
import matplotlib.colors as colors 
import flopy
import matplotlib.pyplot as plt
import numpy as np
from random_fields import import_field
import matplotlib.cm as cm
from matplotlib.colors import LogNorm, LightSource, Normalize
import matplotlib as mpl
import matplotlib.patches as mpatches
from semi_interface import SemiCoastHead
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import utils
from scipy.stats import norm
from random import sample
import scipy.stats as stats

    
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
nstp = 1e6
tsmult = 1
x_onshore = 400
h_onshore = 0.75
anis = 0.01
alpha = 0
kb = 0.1

def plot_3D_field_and_boundary():

    f = plt.figure(figsize=(9.25, 3.5))
    plt.rcParams.update({'font.size': 8})
    gs = f.add_gridspec(2, 6, height_ratios=[5, 0.25])
    ax0 = f.add_subplot(gs[1, :3])
    ax1 = f.add_subplot(gs[1, 3:])
    ax2 = f.add_subplot(gs[0, :3], projection='3d')
    ax3 = f.add_subplot(gs[0, 3:], projection='3d')
    axs = [ax0, ax1, ax2, ax3]

    H = 20
    D = 20
    kb = 0.1
    x_onshore = 400
    x_offshore = 4000
    z_slice = 30
    x_slice = 2000
    width = 1000
    y_slice = 500
    field = import_field("TSim_Out5", H, 10, D, kb, W=1000)
    ls = LightSource(270, 35)
    field = field[0][:int((H + D) / 10), :int(width / 12.5), :]
    norm = LogNorm(vmin=0.1, vmax=np.max(field))
    # plot left base
    x = np.linspace(-x_onshore, x_offshore, int((x_onshore + x_offshore) / 12.5) + 1)
    z = np.linspace(-z_slice, -H - D, int((H + D - z_slice) / 10) + 1)
    X, Z = np.meshgrid(x, z)
    k = field[int(z_slice / 10):, 0, :]
    colors = cm.coolwarm(norm(k))
    ax2.plot_surface(X, 0, Z, facecolors=colors, lightsource=ls)

    # plot left upper
    x = np.linspace(-x_onshore, x_slice, int((x_onshore + x_slice) / 12.5) + 1)
    z = np.linspace(0, -z_slice, int(z_slice / 10) + 1)
    X, Z = np.meshgrid(x, z)
    k = field[:int(z_slice / 10), 0, :int((x_onshore + x_slice) / 12.5)]
    colors = cm.coolwarm(norm(k))
    ax2.plot_surface(X, 0, Z, facecolors=colors, lightsource=ls)

    # plot top left
    x = np.linspace(-x_onshore, x_slice, int((x_onshore + x_slice) / 12.5) + 1)
    y = np.linspace(0, width, int(width / 12.5) + 1)
    X, Y = np.meshgrid(x, y)
    k = field[0, :, :int((x_onshore + x_slice) / 12.5)]
    colors = cm.coolwarm(norm(k))
    Z = np.zeros_like(Y)
    ax2.plot_surface(X, Y, Z, facecolors=colors, lightsource=ls)

    # plot top right
    x = np.linspace(x_slice, x_offshore, int((x_offshore - x_slice) / 12.5) + 1)
    y = np.linspace(y_slice, width, int((width - y_slice) / 12.5) + 1)
    X, Y = np.meshgrid(x, y)
    k = field[0, int(y_slice / 12.5):, int((x_onshore + x_slice) / 12.5):]
    colors = cm.coolwarm(norm(k))
    Z = np.zeros_like(Y)
    ax2.plot_surface(X, Y, Z, facecolors=colors, lightsource=ls)

    # plot middle left
    y = np.linspace(0, y_slice, int(y_slice / 12.5) + 1)
    z = np.linspace(0, -z_slice, int(z_slice / 10) + 1)
    Y, Z = np.meshgrid(y, z)
    k = field[:int(z_slice / 10), :int(y_slice / 12.5), int((x_onshore + x_slice) / 12.5) - 2]
    colors = cm.coolwarm(norm(k))
    ax2.plot_surface(x_slice, Y, Z, facecolors=colors, lightsource=ls)

    # plot middle right
    x = np.linspace(x_slice, x_offshore, int((x_offshore - x_slice) / 12.5) + 1)
    z = np.linspace(0, -z_slice, int(z_slice / 10) + 1)
    X, Z = np.meshgrid(x, z)
    k = field[:int(z_slice / 10), int(y_slice / 12.5), int((x_onshore + x_slice) / 12.5):]
    colors = cm.coolwarm(norm(k))
    ax2.plot_surface(X, y_slice, Z, facecolors=colors, lightsource=ls)

    # plot middle bottom
    x = np.linspace(x_slice, x_offshore, int((x_offshore - x_slice) / 12.5) + 1)
    y = np.linspace(0, y_slice, int(y_slice / 12.5) + 1)
    X, Y = np.meshgrid(x, y)
    k = field[int(z_slice / 12.5), :int(y_slice / 12.5), int((x_onshore + x_slice) / 12.5):]
    colors = cm.coolwarm(norm(k))
    Z = -(z_slice) * np.ones_like(Y)
    ax2.plot_surface(X, Y, Z, facecolors=colors, lightsource=ls)

    # plot right base
    y = np.linspace(0, width, int(width / 12.5) + 1)
    z = np.linspace(-z_slice, -H - D, int((H + D - z_slice) / 10) + 1)
    Y, Z = np.meshgrid(y, z)
    k = field[int(z_slice / 10):, :, -1]
    colors = cm.coolwarm(norm(k))
    ax2.plot_surface(x_offshore, Y, Z, facecolors=colors, lightsource=ls)

    # plot right upper
    y = np.linspace(y_slice, width, int((width - y_slice) / 12.5) + 1)
    z = np.linspace(0, -z_slice, int(z_slice / 10) + 1)
    Y, Z = np.meshgrid(y, z)
    k = field[:int(z_slice / 10), int(y_slice / 12.5):, -1]
    colors = cm.coolwarm(norm(k))
    ax2.plot_surface(x_offshore, Y, Z, facecolors=colors, lightsource=ls)

    ax2.set_xlabel("Distance offshore [m]", labelpad=10.0)
    ax2.set_ylabel("Distance alongshore [m]", labelpad=5.0)
    ax2.set_zlabel("Depth [m]")
    ax2.set_box_aspect((1, 0.8, 0.2))
    ax2.set_zticks([0, -20, -40])
    ax2.set_yticks([0, 500, 1000])
    ax2.set_xticks([0, 1000, 2000, 3000, 4000])
    facies_values = np.sort(np.unique(field), axis=None)
    facies0 = mpatches.Patch(color=cm.coolwarm(norm(facies_values[3])), label='facies 1')
    facies1 = mpatches.Patch(color=cm.coolwarm(norm(facies_values[2])), label='facies 2')
    facies2 = mpatches.Patch(color=cm.coolwarm(norm(facies_values[1])), label='facies 3')
    aquitard = mpatches.Patch(color=cm.coolwarm(norm(facies_values[0])), label='aquitard')
    axs[0].legend(handles=[facies0, facies1, facies2, aquitard], ncol=4, loc="upper center")

    x = np.array([-400, 4000])
    y = np.array([0, 1000])
    X, Y = np.meshgrid(x, y)

    axs[3].plot_wireframe(X, Y, 0*np.ones_like(X), colors='black', lw=0.5)
    axs[3].plot_wireframe(X, Y, -20 * np.ones_like(X), colors='black', lw=0.5)
    axs[3].plot_wireframe(X, Y, -40 * np.ones_like(X), colors='black', lw=0.5)

    z = np.array([0, -40])
    X, Z = np.meshgrid(x, z)

    axs[3].plot_wireframe(X, 0*np.ones_like(X), Z, colors='black', lw=0.5)
    axs[3].plot_wireframe(X, 1000*np.ones_like(X), Z, colors='black', lw=0.5)

    # boundarys
    Y, Z = np.meshgrid(y, z)
    offshore = axs[3].plot_surface(4000, Y, Z, color=cm.viridis(1.0), alpha=0.5, lightsource=ls, label="offshore boundary")
    Y, Z = np.meshgrid(y, np.array([-20, -40]))
    onshore = axs[3].plot_surface(-400, Y, Z, color=cm.viridis(0.0), alpha=0.5, lightsource=ls, label="onshore boundary")
    X, Y = np.meshgrid(np.array([0, 4000]), y)
    axs[3].plot_surface(X, Y, 0*np.ones_like(X), color=cm.viridis(1.0), alpha=0.5, lightsource=ls)
    axs[3].text(1000, 0, -17.5, "Aquitard", "x", va="center")
    axs[3].text(1000, 0, -37.5, "Aquifer", "x", va="center")

    axs[1].legend(handles=[onshore, offshore], ncol=2, loc='upper center')

    axs[3].set_xlabel("Distance offshore [m]", labelpad=10.0)
    axs[3].set_ylabel("Distance alongshore [m]", labelpad=5.0)
    axs[3].set_zlabel("Depth [m]")
    axs[3].set_box_aspect((1, 0.8, 0.2))
    axs[3].set_zticks([0, -20, -40])
    axs[3].set_yticks([0, 500, 1000])
    axs[3].set_xticks([0, 1000, 2000, 3000, 4000])

    axs[0].axis("off")
    axs[1].axis("off")

    f.set_constrained_layout(True)
    f.savefig(f'/home/superuser/objective_2/results/model_setup.png', dpi=600)


def plot_3D_plume(name, i, results=False):
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

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(7, 4))
    ls = LightSource(0, 35)

    
    if results==False:
        name = f"{name}{i}"
        conc, head, qx, qy, qz = modelling.get_last_results(name)   
    elif results == 'steady': 
        results = np.load(f'//collated_outputs/{name}_collated.npy')
        index = np.load(f'//collated_outputs/{name}_collated_index.npy')
        name = f"{name}{i}"
        results = results[index==i]
        conc = results[0]
    elif results == 'transient':
        results = np.load(f'//collated_outputs/{name}_transient_collated.npy')
        index = np.load(f'//collated_outputs/{name}_transient_collated_index.npy')
        name = f"{name}{i}_"
        results = results[index==i]
        conc = results[0]
                  
    plt.rcParams.update({'font.size': 9})

    x = np.linspace(-400, 4000, 352)
    y = np.linspace(0, 1000, 80)
    z = np.linspace(0, -40, 16)
    levels = [0.35, 3.5, 17.75]
    X, Z = np.meshgrid(x, z)
    slices = 3

    viridis = mpl.colormaps['viridis']
    for i in range(slices):
        m = ax.plot_surface(X, i*W/(slices-1), Z, lightsource=ls, facecolors=viridis(conc[:,int(np.floor(i*79/(slices-1))), :]))
        
    ax.set_xlabel("Distance offshore [m]", labelpad=10.0)
    ax.set_ylabel("Distance alongshore [m]", labelpad=5.0)
    ax.set_zlabel("Depth [m]")
    ax.set_box_aspect((1, 0.8, 0.2))
    ax.set_zticks([0, -20, -40])
    ax.set_yticks([0, 500, 1000])   
    ax.set_xticks([0, 1000, 2000, 3000, 4000])
    ax.set_xlim(-400, 4000)

    norm = Normalize(vmin=0, vmax=35)
    m = cm.ScalarMappable(cmap=cm.viridis, norm=norm)
    plt.colorbar(m, ax=ax, location="left", label=r"C [PSU]", shrink=0.75)
    plt.savefig(f'/home/superuser/objective_2/collated_outputs/{name}_collated.png', dpi=600)

def plot_3D_plume_before_and_after(name):
    pre_conc, _, _, _, _ = modelling.get_last_results(name)
    post_conc, _, _, _, _ = modelling.get_last_results(f"{name}_")

    f = plt.figure(figsize=(9, 3.5))
    plt.rcParams.update({'font.size': 8})
    gs = f.add_gridspec(2, 6, height_ratios=[5, 0.25])
    ax0 = f.add_subplot(gs[1, 2:4])
    ax1 = f.add_subplot(gs[0, :3], projection='3d')
    ax2 = f.add_subplot(gs[0, 3:], projection='3d')
    axs = [ax0, ax1, ax2]
    #
    ls = LightSource(0, 35)

    x = np.linspace(-400, 4000, 352)
    y = np.linspace(0, 1000, 80)
    z = np.linspace(0, -40, 16)
    levels = [0.35, 3.5, 17.75]
    X, Z = np.meshgrid(x, z)
    slices = 3

    viridis = mpl.colormaps['viridis']
    for i in range(slices):
        m = axs[1].plot_surface(X, i * W / (slices - 1), Z, lightsource=ls,
                            facecolors=viridis(pre_conc[:, int(np.floor(i * 79 / (slices - 1))), :]))

    viridis = mpl.colormaps['viridis']
    for i in range(slices):
        m1 = axs[2].plot_surface(X, i * W / (slices - 1), Z, lightsource=ls,
                            facecolors=viridis(post_conc[:, int(np.floor(i * 79 / (slices - 1))), :]))

    for ax in axs[1:]:
        ax.set_xlabel("Distance offshore [m]", labelpad=10.0)
        ax.set_ylabel("Distance alongshore [m]", labelpad=5.0)
        ax.set_zlabel("Depth [m]")
        ax.set_box_aspect((1, 0.8, 0.2))
        ax.set_zticks([0, -20, -40])
        ax.set_yticks([0, 500, 1000])
        ax.set_xticks([0, 1000, 2000, 3000, 4000])
        ax.set_xlim(-400, 4000)

    norm = Normalize(vmin=0, vmax=35)
    m = cm.ScalarMappable(cmap=cm.viridis, norm=norm)
    f.colorbar(m, cax=axs[0], location='bottom', label=r"C [PSU]", shrink=0.25)
    f.set_constrained_layout(True)
    # f.show()
    f.savefig(f'/home/superuser/objective_2/results/salinity_before_after_{name}.png', dpi=600)

def plot_all_plumes(name):
    finished_index = np.load(f'collated_outputs/{name}_transient_collated_index.npy')
    for idx, real in enumerate([f"{name}{i}" for i in finished_index]):
        # plot_3D_plume(f"{name}", finished_index[idx], results='transient')
        plot_3D_plume_before_and_after(real)

def plot_one_model(name, results=None):

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

    
    if results==None:
        conc, head, qx, qy, qz = modelling.get_last_results(name)   
        for row in range(80):    
            f, axs = plt.subplots(1, 2, figsize=(4, 6), layout="constrained", sharex=True, sharey=True)
               
            x = np.linspace(-x_onshore, -x_onshore+ncol*delc,ncol)/1000
            y = np.linspace(-D-H, 0, nlay)
            plt.rcParams.update({'font.size': 9})

            for i in range(nlay):
                for j in range(ncol):
                    if abs(conc[i, row, j]) == 1.e30:
                        conc[i, row, j] = np.nan
                        head[i, row, j] = np.nan

            hf_pre = proc.hf(conc[:, row,:], head[:, row,:], nlay, 0, delv)
            qx_pre = qx[::-10, row,::20] / np.sqrt(qx[::-10, row,::20]**2+np.abs(qz[::-10, row,::20])**2)
            qz_pre = qz[::-10, row,::20] / np.sqrt(qx[::-10, row,::20]**2+np.abs(qz[::-10, row,::20])**2)
            axs[0].pcolormesh(x, y, conc[::-1, row,:], cmap="viridis", vmax=35, vmin=0)
            axs[0].quiver(x[::20], y[::10], qx_pre/2, qz_pre/2, color="red", width=0.002, scale=20)
            axs[0].set_box_aspect(0.3)

            head_min = np.nanmin(hf_pre)
            head_max = np.nanmax(hf_pre)
            cp0 = axs[-1].contour(x, y, hf_pre[::-1,:], vmin=head_min, vmax=head_max, colors="black", linewidths=0.75)
            axs[-1].set_box_aspect(0.3)
            axs[0].set_ylabel("Depth [masl]")
            axs[-1].clabel(cp0, inline=True, fontsize=7)
            f.set_constrained_layout(True)
            plt.show()
            # plt.savefig(f"/home/superuser/objective_2/results/{name}_row{row}_final_step.png", dpi=600)

def plot_2D_model_steady(name):

    conc, head, qx, qy, qz = modelling.get_last_results(name)
    f, axs = plt.subplots(1, 2, figsize=(4, 6), layout="constrained", sharex=True, sharey=True)

    x = np.linspace(-x_onshore, -x_onshore + ncol * delc, ncol) / 1000
    y = np.linspace(-D - H, 0, nlay)
    plt.rcParams.update({'font.size': 9})

    for i in range(nlay):
        for j in range(ncol):
            if abs(conc[i, 0, j]) == 1.e30:
                conc[i, 0, j] = np.nan
                head[i, 0, j] = np.nan

    row = 0
    hf_pre = proc.hf(conc[:, row, :], head[:, row, :], nlay, 0, delv)
    qx_pre = qx[::-10, row, ::20] / np.sqrt(qx[::-10, row, ::20] ** 2 + np.abs(qz[::-10, row, ::20]) ** 2)
    qz_pre = qz[::-10, row, ::20] / np.sqrt(qx[::-10, row, ::20] ** 2 + np.abs(qz[::-10, row, ::20]) ** 2)
    axs[0].pcolormesh(x, y, conc[::-1, row, :], cmap="viridis", vmax=35, vmin=0)
    axs[0].quiver(x[::20], y[::10], qx_pre / 2, qz_pre / 2, color="red", width=0.002, scale=20)
    axs[0].set_box_aspect(0.3)
    head_min = np.nanmin(hf_pre)
    head_max = np.nanmax(hf_pre)
    cp0 = axs[-1].contour(x, y, hf_pre[::-1, :], vmin=head_min, vmax=head_max, colors="black", linewidths=0.75)
    axs[-1].set_box_aspect(0.3)
    axs[0].set_ylabel("Depth [masl]")
    axs[-1].clabel(cp0, inline=True, fontsize=7)
    f.set_constrained_layout(True)
    plt.savefig(f"/home/superuser/objective_2/results/{name}_2D_final_step.png", dpi=600)


def results_0(name, n):
    raise NotImplementedError

# plot plan view for each of the models
def results_1(name, n):
   
    interfaces = []
    # for each realization
    # for real in [f"{name}{i}" for i in range(n)]:
    #     # load results
    #     try:
    #         conc, _, _, _, _ = modelling.get_last_results(real) 
    #     # get toe position across the rows
    #         interfaces.append([-x_onshore + delc*np.min(np.argmax([conc[:, row, :] >= 0.35], axis=2)) for row in range(nrow)])
    #     except:
    #         pass
    #     # get toe array
    results = np.load(f'//collated_outputs/{name}_collated.npy')
    for result in results:
        interfaces.append([-x_onshore + delc*np.min(np.argmax([result[:, row, :] >= 0.35], axis=2)) for row in range(nrow)])
    
    f, ax = plt.subplots(figsize=(4, 2))
    ax.plot(np.array(interfaces).T, np.linspace(0, 1000, 80), color="grey")
    ax.set_ylabel("Distance alongshore [m]")
    ax.set_xlabel("Distance offshore [m]")
    f.tight_layout()
    plt.show()

    # plot all
    # raise NotImplementedError

# find the toe position after the head drop    
def results_2(name, n):
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
    nstp = 1e6
    tsmult = 1
    x_onshore = 400
    h_onshore = 0.75
    anis = 0.01
    alpha = 0
    kb = 0.1

    toes = []
    hks = []
    # for each realization
    for i, real in enumerate([f"{name}{i}" for i in range(n)]):
        
        try:
            # load results
            conc, _, _, _, _ = modelling.get_last_results(real) 
            # find toe position 
            toes.append(np.min([-x_onshore + delc*np.min(np.argmax([conc[:, row, :] >= 0.35], axis=2)) for row in range(nrow)]))
            
            missing = False
        except:
            missing = True

        if not missing:
            hk_eff, hv_eff = modelling.get_bulk_conductivities(L, H, W, i)
            hks.append(hk_eff)
    
    f, ax = plt.subplots(figsize = (3,3)) 
    ax.scatter(hks, toes, c="grey") 
    ax.set_xlabel(r"$K_a^{eff}$ [m/day]")
    ax.set_ylabel("Distance offshore [m]")
    ax.set_box_aspect(1)
    f.tight_layout()
        # get toe arra
    # for each realization
        # load result
        # find toe position 
        # find bulk hydraulic conductivity
    # plot all
    plt.show()
    #   

# find maximum onshore extent and after head drop
def results_3(name, n):
       
    D = 20
    H = 20
    W = 100
    L = 4400
    delc = 6.25
    delv = 2.25
    delr = 6.25
    nlay = int((D+H)/delv)
    ncol = int(L/delc)
    nrow = int(W/delr)
    perlen = 1e8
    nper = 1
    nstp = 1e6
    tsmult = 1
    x_onshore = 400
    h_onshore = 0.75
    anis = 0.01
    alpha = 0
    kb = 0.1
    pre_toes = []
    post_toes = []

    for real in [f"{name}{i}" for i in range(n)]:
        # load results
        try:
            ws = f"./model_files/{real}"
            sim = flopy.mf6.MFSimulation.load(
                sim_ws=ws, exe_name= "/home/superuser/mf6", verbosity_level=0,)
            strt_head = list(sim._models.values())[0].ic.strt.data
            strt_conc = list(sim._models.values())[1].ic.strt.data
            conc, _, _, _, _ = modelling.get_last_results(real) 

            pre_toes.append(np.min([-x_onshore + delc*np.min(np.argmax([strt_conc[:, row, :] >= 0.35], axis=2)) for row in range(nrow)]))
            post_toes.append(np.min([-x_onshore + delc*np.min(np.argmax([conc[:, row, :] >= 0.35], axis=2)) for row in range(nrow)]))


        except:
            pass
    # for each realization
        # load result
        # find toe position 

    f, ax = plt.subplots(figsize=(3,3))
    ax.set_box_aspect(1)
    ax.scatter(pre_toes, post_toes)
    ax.plot(np.linspace(np.min(pre_toes)-50, np.max(pre_toes)+50), np.linspace(np.min(pre_toes)-50, np.max(pre_toes)+50), color="grey", ls="--")
    ax.set_xlabel("Distance offshore [m]")
    ax.set_ylabel("Distance offshore [m]")
    f.tight_layout()
    plt.show()
    # plot all
    raise NotImplementedError

def steady_results_multiplot(name):

    W = 1000
    delc = 12.5
    delr = 12.5
    nrow = int(W/delr)

    x_onshore = 400

    f = plt.figure(figsize=(8,5))
    plt.rcParams.update({'font.size': 9})
    gs = f.add_gridspec(3, 3, height_ratios=[5, 5, 1])
    ax0 = f.add_subplot(gs[0, 1:])
    ax1 = f.add_subplot(gs[1, 0])
    ax2 = f.add_subplot(gs[1, 1])
    ax3 = f.add_subplot(gs[1, 2])
    ax_leg0 = f.add_subplot(gs[0, 0])
    ax_bottom_leg = f.add_subplot(gs[2, :])
    twin = ax_leg0.twinx()

    interfaces = []
    results = np.load(f'collated_outputs/{name}_collated.npy')
    for result in results:
        interfaces.append([-x_onshore + delc*np.min(np.argmax([result[:, row, :] >= 0.35], axis=2)) for row in range(nrow)])
    
    ax0.plot(np.array(interfaces)[0].T, np.linspace(0, 1000, 80), color="grey", lw=0.5, label='realizations')
    ax0.plot(np.array(interfaces)[1:].T, np.linspace(0, 1000, 80), color="grey", lw=0.5, label=None)

    averages = np.average(np.array(interfaces).T, 0)
    average_idx = np.argwhere(averages == np.median(averages))[0][0]
    ax0.plot(interfaces[average_idx], np.linspace(0, 1000, 80), color='red', label='median')
    ax0.set_ylabel("distance alongshore [m]")
    ax0.set_xlabel(r"$x_T$ [m]")
    h, l = ax0.get_legend_handles_labels()
    ax_leg0.legend(h, l, loc='upper left', title='Legend (right)')

    widths = []
    for result in results:
        widths.append(-delc*(np.min(np.argmax(result[:, :, :] >= 0.35, axis=2)-
                                   np.max(np.argmax(result[:, :, :] >= 0.35, axis=2)))))

    ax2.hist(np.average(np.array(interfaces), axis=1), bins=12, color='grey', density=True, alpha=0.5)
    ax2.set_ylabel('density')
    ax2.set_xlabel(r'$x_T^{mean}$ [m]')

    hks = np.load(f'collated_outputs/{name}_hk_eff.npy')
    hks_analytical = np.linspace(np.min(hks), np.max(hks), 100)
    xTs_analytical  = []
    for hk in hks_analytical:
        xTs_analytical.append(SemiCoastHead(k=hk, H=20, c=20/0.001, h=0.75, x=-400, 
                                            rhof=1000, rhos=1025, Ls=4000, ztop=-20, sealevel=0).xtoe)

    ax1.plot(hks_analytical, xTs_analytical, color='blue', ls='--', label=r"$x_T$ (Bakker, 2017)")
    ax1.scatter(hks, np.average(np.array(interfaces), axis=1), marker='+', linewidths=0.75, c='grey', s=30, label=r'$x_T^{mean}$', zorder=4)
    ax1.set_xlabel(r'$K_a^{eff}$ [m/day]')
    ax1.set_ylabel('distance offshore [m]')
    h, l = ax1.get_legend_handles_labels()
    twin.legend(h, l, loc="lower left", title='Legend (below)')
    
    twin.axis("off")
    ax_leg0.axis("off")


    average_widths = []
    for result in results:
        extent = np.argmax(result[:, :, :] >= 31.5, axis=2)
        extent[extent == 0] = ncol
        average_widths.append(delc*np.mean(extent - np.argmax(result[:, :, :] >= 3.5, axis=2)))

    ax3.hist(average_widths, bins=10, density=True, alpha=0.5, color='grey', label="realizations")
    ax3.set_ylabel('density')
    ax3.set_xlabel(r'mean interface width [m]')

    steady_stats = utils.find_steady_description(name)
    norm_toe = norm(steady_stats[0], steady_stats[1])
    norm_width = norm(steady_stats[2], steady_stats[3])

    x = np.linspace(0, 3000, 3000)
    ax2.plot(norm_toe.pdf(x), color="orange")
    ax3.plot(norm_width.pdf(x), color="orange", label="estimated pdf")

    h,l = ax3.get_legend_handles_labels()
    ax_bottom_leg.legend(h, l, loc='upper right', ncol=2)

    ax_bottom_leg.axis("off")
    f.set_constrained_layout(True)
    # plt.show()
    f.savefig(f'/home/superuser/objective_2/results/{name}_steady_multiplot.png', dpi=600)
    
def plot_interface_movement_plan_view(name):

    pre = np.load(f'//collated_outputs/{name}_collated.npy')
    post = np.load(f'//collated_outputs/{name}_transient_collated.npy')
    pre_index = np.load(f'//collated_outputs/{name}_collated_index.npy')
    post_index = np.load(f'//collated_outputs/{name}_transient_collated_index.npy')

    for realization_post, index in zip(post, post_index):
        realization_pre = pre[np.argwhere(pre_index == index)]
        f, ax = plt.subplots(1,1, figsize=(6, 3))
        interface_post = [-x_onshore + delc*np.min(np.argmax([realization_post[:, row, :] >= 0.35], axis=2)) for row in range(nrow)]
        interface_pre = [-x_onshore + delc*np.min(np.argmax([realization_pre[0][0][:, row, :] >= 0.35], axis=2)) for row in range(nrow)]
        ax.plot(np.array(interface_pre).T, np.linspace(0, 1000, 80), ls='--', color="grey", label=r'$x_T^{pre}$')
        ax.plot(np.array(interface_post).T, np.linspace(0, 1000, 80), color="grey", label=r'$x_T^{post}$')
        ax.fill_betweenx(np.linspace(0, 1000, 80), interface_pre, interface_post, color='grey', alpha=0.5, label='area salinized')
        ax.set_xlim(-400, 3000)
        ax.set_ylabel('distance alongshore [m]')
        ax.set_xlabel('distance offshore [m]')
        ax.set_title("t = 100 years, h = -1 m")
        ax.legend()
        f.savefig(f'/home/superuser/objective_2/results/interface_movement{index}.png', dpi=600)

def transient_results_multi(name):

    f = plt.figure(figsize=(8,5))
    plt.rcParams.update({'font.size': 9})
    gs = f.add_gridspec(3, 3, height_ratios=[5, 5, 1])
    ax0 = f.add_subplot(gs[0, 1:])
    ax1 = f.add_subplot(gs[1, 0])
    ax2 = f.add_subplot(gs[1, 1])
    ax3 = f.add_subplot(gs[1, 2])
    ax_bottom_leg = f.add_subplot(gs[2, :])
    ax_leg0 = f.add_subplot(gs[0, 0])
    twin = ax_leg0.twinx()

    interfaces = []
    results = np.load(f'collated_outputs/{name}_transient_collated.npy')
    for result in results:
        interfaces.append([-x_onshore + delc*np.min(np.argmax([result[:, row, :] >= 0.35], axis=2)) for row in range(nrow)])
    
    ax0.plot(np.array(interfaces)[0].T, np.linspace(0, 1000, 80), color="grey", lw=0.5, label='realizations')
    ax0.plot(np.array(interfaces)[1:].T, np.linspace(0, 1000, 80), color="grey", lw=0.5, label=None)
    ax0.set_ylabel("Distance alongshore [m]")
    ax0.set_xlabel("Distance offshore [m]")
    averages = np.average(np.array(interfaces).T, 0)
    average_idx = np.argwhere(averages == np.median(averages))[0][0]
    ax0.plot(interfaces[average_idx], np.linspace(0, 1000, 80), color='red', label='median')
    ax0.set_ylabel("distance alongshore [m]")
    ax0.set_xlabel(r"$x_T$ [m]")
    h, l = ax0.get_legend_handles_labels()
    ax_leg0.legend(h, l, loc='upper left', title='Legend (right)')

    max = np.load(f'collated_outputs/{name}_transient_collated_max_toe.npy')
    average = np.load(f'collated_outputs/{name}_transient_collated_average_toe.npy')
    width = np.load(f'collated_outputs/{name}_transient_collated_average_width.npy')

    t = np.linspace(0, 100, 11)
    ax2.fill_between(t, np.percentile(max, 0, axis=0), np.percentile(max, 100, axis=0), color='grey', alpha=0.25)
    ax2.fill_between(t, np.percentile(max, 25, axis=0), np.percentile(max, 75, axis=0), color='grey', alpha=0.25)
    ax2.plot(t, np.median(max, axis=0), color='grey')
    ax2.set_ylabel(r"$x_T^{min}$ [m]")
    ax2.set_xlabel("time [years]")

    t = np.linspace(0, 100, 11)
    ax3.fill_between(t, np.percentile(width, 0, axis=0), np.percentile(width, 100, axis=0), color='grey', alpha=0.25)
    ax3.fill_between(t, np.percentile(width, 25, axis=0), np.percentile(width, 75, axis=0), color='grey', alpha=0.25)
    ax3.plot(t, np.median(width, axis=0), color='grey')
    ax3.set_ylabel("mean interface width [m]")
    ax3.set_xlabel("time [years]")

    index = np.load(f'collated_outputs/{name}_transient_collated_index.npy')
    hks = np.load(f'collated_outputs/{name}_hk_eff.npy')
    hks_analytical = np.linspace(np.min(hks), np.max(hks), 100)
    
    v_T_analytical = []
    porosity = 0.5
    h_f = ((1025-1000)/1000)*(20)
    delta_h = -(-1-h_f)
    for hk in hks_analytical:
        v_T_analytical.append(365*delta_h*hk/(4400*porosity))
        
    difference_average = np.diff(average, axis=1)
    difference_max = np.diff(max, axis=1)
    difference_average[difference_max == 0.0] = np.nan
    difference_max[difference_max== 0.0] = np.nan
    velocity = difference_average/10
    velocity_max = difference_max/10

    #  v_Ts = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_transient_v_T.npy')
    ax1.plot(hks_analytical, v_T_analytical, color='blue', ls='--', label=r"$v_T$ (Werner, 2017)")
    ax1.scatter(hks[index[:-1]], -np.nanmean(velocity_max, axis=1)[:-1], marker='+', linewidths=0.75, c='grey', s=30, label=r'$dx_T^{min}$/$dt$', zorder=4)
    ax1.set_xlabel(r'$K_a^{eff}$ [m/day]')
    ax1.set_ylabel(r'velocity [m/year]')
    h, l = ax1.get_legend_handles_labels()
    twin.legend(h, l, loc="lower left", title='Legend (below)')

    patch_outer = mpatches.Patch(color='grey', label='percentiles 0-100', alpha=0.25)
    patch_inner = mpatches.Patch(color='grey', label='percentiles 25-75', alpha=0.5)
    median = Line2D([0], [0], color="grey", label="median")
    ax_bottom_leg.legend(handles=[patch_outer, patch_inner, median], loc='upper right', ncol=3)

    ax_bottom_leg.axis("off")
    twin.axis("off")
    ax_leg0.axis("off")
    f.set_constrained_layout(True)
    f.savefig(f'/home/superuser/objective_2/results/{name}_transient_multiplot.png', dpi=600)

def plot_plan_view_examples(name):
    pre = np.load(f'collated_outputs/{name}_collated.npy')
    post = np.load(f'collated_outputs/{name}_transient_collated.npy')
    pre_index = np.load(f'collated_outputs/{name}_collated_index.npy')
    post_index = np.load(f'collated_outputs/{name}_transient_collated_index.npy')

    plot = sample(list(post_index), 9)
    cmap = cm.get_cmap('Set1')
    colors = [cmap(i/9) for i in range(len(plot))]
    f, axs = plt.subplots(4, 3, figsize=(8, 5), gridspec_kw={'hspace': 0, 'wspace': 0})
    plt.rcParams.update({'font.size': 8})
    for i, idx in enumerate(plot):
        realization_pre = pre[np.argwhere(pre_index == idx)][0][0]
        realization_post = post[np.argwhere(post_index == idx)][0][0]

        interface_post = [-x_onshore + delc * np.min(np.argmax([realization_post[:, row, :] >= 0.35], axis=2)) for row
                          in range(nrow)]
        interface_pre = [-x_onshore + delc * np.min(np.argmax([realization_pre[:, row, :] >= 0.35], axis=2)) for
                         row in range(nrow)]

        axs[i//3][i%3].plot(np.array(interface_pre).T, np.linspace(0, 1000, 80), ls='--', color=colors[i], label=r'$x_T^{pre}$')
        axs[i//3][i%3].plot(np.array(interface_post).T, np.linspace(0, 1000, 80), color=colors[i], label=r'$x_T^{post}$')
        axs[i//3][i%3].fill_betweenx(np.linspace(0, 1000, 80), interface_pre, interface_post, color=colors[i], alpha=0.5, label='area salinized')
        axs[i//3][i%3].set_xlim(-400, 3000)
        axs[i//3][i%3].set_xticklabels([])
        axs[i//3][i%3].set_yticklabels([])
        t = axs[i//3][i%3].text(0.95, 0.85, f'realization {idx}', horizontalalignment='right', verticalalignment='center',
                 transform=axs[i//3][i%3].transAxes)
        t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='None'))

    patch = mpatches.Patch(color='grey', label='area salinzed', alpha=0.5)
    final = Line2D([0], [0], color="grey", label="final interface")
    start = Line2D([0], [0], color="grey", ls="--", label="steady interface")
    axs[-1][-1].legend(handles=[start, final, patch], loc='center right')
    axs[-2][1].set_xlabel("offshore position")
    axs[1][0].set_ylabel("alongshore position")

    axs[-1][-1].axis("off")
    axs[-1][0].axis("off")
    axs[-1][1].axis("off")
    # f.set_constrained_layout(True)
    f.savefig(f'/home/superuser/objective_2/results/{name}_movement_realizations.png', dpi=600)
    plt.show()
    pass

def plot_priors(name, priors, units, types):
    """

    :param name: name of 2d model (dispersive or dual porosity)
    :param priors:
    :param units:
    :param types: types of each of the prior distributions
    """
    f, axs = plt.subplots(1, len(priors), figsize=(2*len(priors), 1.5))
    for i, (param, dist) in enumerate(priors.items()):

        if types[i] == 'uniform':
            vmin = dist.dist.kwds['loc'] - 1
            vmax = dist.dist.kwds['loc'] + dist.dist.kwds['scale'] + 1
        elif types[i] == 'norm':
            vmin = dist.dist.kwds['loc'] - 3*dist.dist.kwds['scale']
            vmax = dist.dist.kwds['loc'] + 3*dist.dist.kwds['scale']

        sample = np.linspace(vmin, vmax, 1000)
        axs[i].plot(sample, dist.dist.pdf(sample), color="orange")
        axs[i].set_ylabel('probability')
        axs[i].set_xlabel(f'{param}, [log({units[i]})]')

    f.set_constrained_layout(True)
    plt.savefig(f'/home/superuser/objective_2/results/priors_{name}.png', dpi=600)

def plot_posteriors():
    pass

if __name__=="__main__":
    # pass
    plot_3D_field_and_boundary()
    # plot_all_plumes("heta03Dc")
    #steady_results_multiplot("heta03Dc")
    #plot_interface_movement_plan_view("heta03Dc")
    #transient_results_multi("heta03Dc")
    # plot_plan_view_examples("heta03Dc")
