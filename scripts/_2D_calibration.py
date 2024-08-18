from scipy.stats import uniform, norm
from pydream.parameters import SampledParam
from pydream.core import run_dream
import numpy as np
import plots
import scipy.stats as stats
import simple_modelling
import modelling
import tempfile
from pathlib import Path
from komanawa.modeltools.model_tools_package.utils.convergence_check import modflow_converged
import utils
from scripts.plots import plot_posteriors


def log_likelihood_dispersive(params):
    k, alpha = params
    with tempfile.TemporaryDirectory() as tdir:
        tdir = Path(tdir)
        base_params = {'ka': 10**k}
        transport_kwargs = {'alpha': 10**alpha}
        swt = simple_modelling.build_2D_steady_model(tdir, 'temp', 'dispersive', base_params, transport_kwargs)
        modelling.run_model(swt)

        toe, width, conc, head = simple_modelling.get_steady_metrics_and_arrays(tdir, 'temp')

        if modflow_converged(tdir.joinpath('mfsim.lst')):
            results_3D = utils.find_steady_description("heta03Dc")
            deltas = [results_3D[0]-toe, results_3D[2]-width]
            sigmas = [results_3D[1], results_3D[3]]
            likelihood = -len(deltas) / 2 * np.log(2 * np.pi) - (np.log(sigmas)).sum() - 1 / 2 * (
                    (np.array(deltas)/np.array(sigmas)) ** 2).sum()
            np.save(f"temp_results/conc_dispersive_k={k}_alpha={alpha}.npy", conc)
            np.save(f"temp_results/head_dispersive_k={k}_alpha={alpha}.npy", head)
            np.save(f"temp_results/metrics_dispersive_k={k}_alpha={alpha}.npy", [toe, width])
        else:
            likelihood = -np.inf
            # save final head and conc
        return likelihood


def log_likelihood_dual(params):
    k, volfrac, zeta_im = params
    with tempfile.TemporaryDirectory() as tdir:
        tdir = Path(tdir)
        base_params = {'ka': k**2}
        transport_kwargs = {'alpha': 0.01, 'porosity_im': 0.5, 'volfrac': volfrac**2, 'zetaim': zeta_im**2}
        swt = simple_modelling.build_2D_steady_model(tdir,'temp', 'dual porosity', base_params, transport_kwargs)
        modelling.run_model(swt)

        toe, width, conc, head = simple_modelling.get_steady_metrics_and_arrays(tdir, 'temp')

        if modflow_converged(tdir.joinpath('mfsim.lst')):
            results_3D = utils.find_steady_description("heta03Dc")
            deltas = [results_3D[0]-toe, results_3D[2]-width]
            sigmas = [results_3D[1], results_3D[3]]
            likelihood = -len(deltas) / 2 * np.log(2 * np.pi) - (np.log(sigmas)).sum() - 1 / 2 * (
                    (np.array(deltas)/np.array(sigmas)) ** 2).sum()
            np.save(f"temp_results/conc_dual_k={k}_volfrac={volfrac}_zeta_im={zeta_im}.npy", conc)
            np.save(f"temp_results/head_dispersive_k={k}_volfrac={volfrac}_zeta_im={zeta_im}.npy.npy", head)
            np.save(f"temp_results/metrics_dispersive_k={k}_volfrac={volfrac}_zeta_im={zeta_im}.npy.npy", [toe, width])
        else:
            likelihood = -np.inf
            # save final head and conc
        return likelihood

def generate_priors(name, type, plot=False):
    priors = {}
    if type == "dispersive":
        # log hydraulic conductivity
        ka = SampledParam(norm, loc=1.3, scale=0.4)
        # log dispersivity 0.1 --> 100
        alpha = SampledParam(uniform, loc=-1, scale=3)
        priors['ka'] = ka
        priors["alpha"] = alpha
        units = ["m/day", "m"]
        types = ['norm', 'uniform']

    elif type == "dual porosity":
        # more info about this sort of thing in Sarris et al. (2018)
        ka = SampledParam(norm, loc=1.3, scale=0.4)
        volfrac = SampledParam(uniform, loc=-3, scale=2.4)
        zeta_im = SampledParam(uniform, loc=-4, scale=3) # from Bianchi and Zheng (2016)
        priors['ka'] = ka
        priors["volfrac"] = volfrac
        priors["zeta_im"] = zeta_im
        units = ["m/day", "-", "1/day"]
        types = ['norm', 'uniform', 'uniform']

    if plot:
        plots.plot_priors(type, priors, units, types)
        plot_posteriors(name, 0.5, priors, units, types)

    return priors

def calibrate_steady_state(name, type='dispersive'):
    params = generate_priors(type=type)
    model_name = Path('./dreamz_outputs')
    model_name.mkdir(exist_ok=True)
    model_name = model_name.joinpath(name)
    model_name.mkdir(exist_ok=True)
    if type == 'dispersive':
        params = [params["ka"], params["alpha"]]
        sampled_params, log_ps = run_dream(likelihood=log_likelihood_dispersive, parameters=params, nchains=8,
                                           niterations=100, model_name=f"{str(model_name)}/")
        np.save(f"{model_name}/sampled_params.npy", sampled_params)
        np.save(f"{model_name}/log_ps.npy", log_ps)
        # save sampled_params
        # save log_ps
        # collate and move results

    if type == 'dual porosity':
        params = [params["ka"], params["volfrac"], params["zeta_im"]]
        sampled_params, log_ps = run_dream(likelihood=log_likelihood_dual, parameters=params, nchains=8,
                                           niterations=100, model_name=f"{str(model_name)}/")
        np.save(f"{model_name}/sampled_params.npy", sampled_params)
        np.save(f"{model_name}/log_ps.npy", log_ps)


def test_steady_realization(type="dispersive"):

    if type == "dispersive":
        priors = generate_priors("dispersive", plot=False)
        base_params = {'ka': 20}
        transport_kwargs = {'alpha': 10}
        swt = simple_modelling.build_2D_steady_model("tests", 'dispersive_l', 'dispersive', base_params, transport_kwargs)
        modelling.run_model(swt)

    elif type == "dual porosity":
        priors = generate_priors("dual porosity", plot=False)
        base_params = {'ka': 20}
        transport_kwargs = {'alpha': 0.01, 'porosity_im':0.5, 'volfrac':0.1, 'zetaim':0.005}
        swt = simple_modelling.build_2D_steady_model('dual_s', 'dual porosity', base_params, transport_kwargs)
        modelling.run_model(swt)



if __name__ == "__main__":
    generate_priors(name='test_dispersive', type="dispersive", plot=True)
    generate_priors(name='test_dual', type="dual porosity", plot=True)
    # test_steady_realization('dispersive')
    # plots.plot_2D_model_steady('dispersive_l')
    # calibrate_steady_state("test_dispersive", type="dispersive")
    # calibrate_steady_state("test_dual", type="dual porosity")
    # log_likelihood_dispersive([10, 10])