import modelling
import plots
from multiprocessing import Pool
import numpy as np

def test_one_model(inputs):
    name, field_name = inputs
    sim = modelling.build_base_model(name, field_name)
    modelling.run_model(sim)
    
def build_and_run_base_models(name, n):


    ran = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_collated_index.npy')

    inputs = []
    for i in range(n):
        if i not in ran:
            inputs.append([f"{name}{i}", f"TSim_Out{i+1}"])

    p=Pool(processes=16)
    p.map(test_one_model, inputs)
    #test_one_model(inputs[0]) 

def run_reduced_head(inputs):
    name, field_name = inputs
    try:    
        conc, head, _, _, _ = modelling.get_last_results(name)
        sim = modelling.add_head_to_model(f"{name}_", field_name, conc, head)
        modelling.run_model(sim)
    except:
        pass

def build_and_run_modern(name, n):
    
    ran = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_index.npy')

    inputs = []
    for i in range(n):
        if i not in ran:
            inputs.append([f"{name}{i}", f"TSim_Out{i+1}"]) 

    p=Pool(processes=4)
    p.map(run_reduced_head, inputs)

def run_pumping():
    pass

def run_recovery():
    pass

def run_modpath_for_all_modern(name):
    s_index = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_collated_index.npy')
    t_index = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_transient_collated_index.npy')
    interfaces = np.load(f'/home/superuser/objective_2/collated_outputs/{name}_collated_interface_cols.npy')
    runs = []

    for i in t_index:
        real_name = f"{name}{i}_"
        interface = interfaces[np.argwhere(s_index == i)[0][0]]
        runs.append([real_name, interface])

    # p = Pool(processes=4)
    # p.map(modelling.run_modpath_modern, runs)
    modelling.run_modpath_modern(runs[1])

if __name__ == "__main__":
    # build_and_run_base_models("heta03Dc", 30)
    #test_one_model("het1a03d", "TSim_Out1")
    # plots.plot_one_model("heta03Dc1")
    # plots.plot_3D_plume("heta03Dc1")
    # plots.results_1("heta03Dc", 30)
    # plots.results_2("heta03Dc", 10)
    # build_and_run_modern("heta03Dc", 30)
    # plots.results_3("heta03Dfine", 10)
    run_modpath_for_all_modern("heta03Dc")