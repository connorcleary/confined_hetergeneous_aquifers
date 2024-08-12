import flopy
import numpy as np
import modelling
import post_processing as proc

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

def collate_transient_results(name, n):
    results = []
    finished_index = []
    for idx, real in enumerate([f"{name}{i}_" for i in range(n)]):
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
        conc = modelling.get_all_conc(f'{name}{idx}_')
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

if __name__=="__main__":
    # collate_transient_results('heta03Dc', 30)
    # find_effective_conductivities('heta03Dc')
   # find_toe_position_over_time
   # find_toe_position_over_time
    find_metrics_over_time('heta03Dc')
    # find_steady_description('heta03Dc')