import flopy
import numpy as np
import modelling
import post_processing as proc
from pathlib import Path

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

def collate_steady_results(name, n):
    results = []
    finished_index = []
    for idx, real in enumerate([f"{name}{i}" for i in range(n)]):
        conc, _, _, _, _ = modelling.get_last_results(real) 
        if conc is not None:
            finished_index.append(idx)
            results.append(conc)

    np.save(f'/home/superuser/objective_2/collated_outputs/{name}_collated.npy', results)
    np.save(f'/home/superuser/objective_2/collated_outputs/{name}_collated_index.npy', finished_index)


def collate_steady_interface_cells(name):
    concs = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_collated.npy')
    interfaces = []
    for idx, conc in enumerate(concs):
        interfaces.append([np.min(np.argmax([conc[:, row, :] >= 0.35], axis=2)) for row in range(conc.shape[1])])

    np.save(f'/home/superuser/objective_2/collated_outputs/{name}_collated_interface_cols.npy', interfaces)


def collate_transient_results(name, n):
    results = []
    finished_index = []
    for idx, real in enumerate([f"{name}{i}_t" for i in range(n)]):
        conc, _, _, _, _ = modelling.get_last_results(real) 
        if conc is not None:
            finished_index.append(idx)
            results.append(conc)

    np.save(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated.npy', results)
    np.save(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_index.npy', finished_index)

def find_metrics_over_time(name):
    finished_index = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_index.npy')
    max_toe = []
    average_toe = []
    average_width = []

    for idx in finished_index:

        conc0 = modelling.get_last_results(f'{name}{idx}')[0]
        conc = modelling.get_all_conc(f'{name}{idx}_t')
        temp_max = [proc.find_max_toe_position(conc0, delc, x_onshore)]
        temp_ave = [proc.find_average_toe_position(conc0, delc, x_onshore)]
        temp_width = [proc.find_average_width(conc0, delc, ncol, x_onshore)]

        for conc_t in conc:
            temp_max.append(proc.find_max_toe_position(conc_t, delc, x_onshore))
            temp_ave.append(proc.find_average_toe_position(conc_t, delc, x_onshore))
            temp_width.append(proc.find_average_width(conc_t, delc, ncol, x_onshore))

        max_toe.append(temp_max)
        average_toe.append(temp_ave)
        average_width.append(temp_width)

    np.save(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_max_toe.npy', max_toe)
    np.save(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_average_toe.npy', average_toe)
    np.save(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_average_width.npy', average_width)

def find_velocity_over_time(name):
    finished_index = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_index.npy')
    velocity = []
    t = [i * 3600 / 365 for i in range(11)] + [100]

    for idx in finished_index:

        conc0 = modelling.get_last_results(f'{name}{idx}')[0]
        conc = modelling.get_all_conc(f'{name}{idx}_t')
        interfaces = [[-x_onshore + delc*np.min(np.argmax([conc0[:, row, :] >= 0.35], axis=2)) for row in range(conc0.shape[1])]]
        for conc_t in conc:
            interfaces.append([-x_onshore + delc*np.min(np.argmax([conc_t[:, row, :] >= 0.35], axis=2)) for row in range(conc_t.shape[1])])
        interfaces = np.asarray(interfaces)
        diff = np.diff(interfaces, axis=0)
        velocity.append(np.divide(diff.T, np.diff(t)))

    np.save(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_velocity.npy', velocity)


def collate_shoreline_salinity(name):
    finished_index = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_index.npy')
    results = []
    for idx in finished_index:
        conc0 = modelling.get_last_results(f'{name}{idx}')[0]
        conc = modelling.get_all_conc(f'{name}{idx}_t')
        conc = np.asarray(conc)
        all_conc = np.concatenate(([conc0], conc), axis=0)[:, :, :, int(400/12.5)]
        all_conc = np.max(all_conc, axis=1)
        results.append(all_conc)

    np.save(f'/home/superuser/objective_2/collated_outputs/{name}_average_shoreline_conc.npy', results)


def find_effective_conductivities(name):

    H = 20
    W = 1000
    L = 4400
    index = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_collated_index.npy')

    hks = []
    hvs = []
    # for each realization
    for idx in index:
        hk_eff, hv_eff = modelling.get_bulk_conductivities(L, H, W, idx)
        hks.append(hk_eff)
        hvs.append(hv_eff)

    np.save(f'/home/superuser/objective_2/collated_outputs/{name}_hk_eff.npy', hks)
    np.save(f'/home/superuser/objective_2/collated_outputs/{name}_hv_eff.npy', hvs)

def find_steady_description(name):
    average_toe = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_average_toe.npy')[:,0]
    average_width = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_average_width.npy')[:,0]
    toe_mean = np.average(average_toe)
    toe_std = np.sqrt(np.var(average_toe))
    width_mean = np.average(average_width)
    width_std = np.sqrt(np.var(average_width))
    return [toe_mean, toe_std, width_mean, width_std]

def find_metrics_over_time_2D():
    pass

def collate_pathline_data(name):
    index = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_index.npy')
    results = []
    # aquifer only
    # one particle in upper and lower layer
    # 20 particles in total
    rows = np.linspace(0, 79, 40)

    particles = []
    for i, row in enumerate(rows):
        if i%2 == 0:
            particles.append(int(row*8+2))
        if i%2 == 1:
            particles.append(int(row*8+6))

    for idx in index:
        result = []
        file = f"/home/superuser/objective_2/model_files/{name}{idx}_t/{name}{idx}_t_mp.mppth"
        try:
            p = flopy.utils.PathlineFile(file)
            rec_array = p.get_alldata()
            for particle in particles:
                ts = []
                for step in rec_array[particle]:
                    ts.append([step[0], step[1]])
                result.append(ts)
            results.append(result)
        except:
            results.append(None)

    return results

def get_n_best(name, nplot):
    """
    Function to get the n best runs from a dreamz output
    :param name:
    :param nplot:
    :return: sampled_params: all the sampled parameters
    :return: best: the indices of the best n runs
    """
    # get the index of the n best runs
    dreamz_dir = Path(f'/home/superuser/objective_2/dreamz_outputs/{name}/')
    sampled_params = np.load(dreamz_dir.joinpath('sampled_params.npy'))
    log_ps = np.load(dreamz_dir.joinpath('log_ps.npy'))
    log_ps = log_ps.flatten()
    sampled_params = sampled_params.reshape(len(log_ps), sampled_params.shape[-1])
    best = np.argsort(log_ps)[-int(nplot*len(log_ps)):]

    return sampled_params[best], best


if __name__=="__main__":
    pass
    # collate_transient_results('heta03Dc', 30)
    # find_effective_conductivities('heta03Dc')
   # find_toe_position_over_time
   # find_toe_position_over_time
    # find_metrics_over_time('heta03Dc')
    # find_steady_description('heta03Dc')
    # collate_shoreline_salinity('heta03Dc')
    # collate_steady_interface_cells('heta03Dc')
    #collate_transient_results('heta03Dc', 30)
    # collate_pathline_data('heta03Dc')
    find_velocity_over_time('heta03Dc')