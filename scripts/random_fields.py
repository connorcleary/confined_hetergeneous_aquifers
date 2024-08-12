import numpy as np

def import_field(field_name, H, delv, D, kb, W=3600, delr=12.5, L=4400, delc=12.5):
    # import, change resolution, change to facies and add aquitard
    with open(f'./TW_3a_cond_hydrofacies_realisations/hydrofacies_categories/{field_name}.grd', 'r') as file:
        ncol, nrow, nlay = map(int, file.readline().split())

        # Read the entire file into a flattened NumPy array
        iarray_flat = np.fromfile(file, dtype=float, count=-1, sep=' ')

        # Reshape the flattened array into a 3D array. Fortran-style ordering (column-major)
        iarray = iarray_flat.reshape((ncol, nrow, nlay), order='F')
        # reshape for flopy
        iarray = np.abs(iarray.transpose((2, 1, 0)))
        # truncate
        iarray = iarray[:int((H/100)*nlay),:,:]
        iarray = iarray[:,:int((W/3600)*nrow),:]
        iarray = iarray[:,:,:int((L/4400)*ncol)]

        # import hydrofacies properties
        i = int(field_name.strip("TSim_Out")) - 1
        facies_k = np.genfromtxt('./TW_3a_cond_hydrofacies_realisations/hydrofacies_conductivity_list.csv', skip_header=1, delimiter=",")[i, 2:]
        # change from facies to k
        iarray[iarray==1] = facies_k[0] 
        iarray[iarray==2] = facies_k[1]
        iarray[iarray==3] = facies_k[2]
        # increase_frequency of sample
        iarray = iarray.repeat(int(10/delv), axis=0)
        iarray = iarray.repeat(int(12.5/delr), axis=1)
        iarray = iarray.repeat(int(12.5/delc), axis=2)
        # add aquitard
        hk = np.concatenate((kb*np.ones((int(D/delv), int(W/delr), int(L/delc))), iarray), axis=0)
        hka = [0.01 for i in range(int(D/delv))] +  [1 for i in range(int(H/delv))] #  # [13/50 for i in range(int(H/delv))]#
        return hk, hka

def import_bulk_field(field_name, H, delv, D, kb, W=3600, delr=12.5, L=4400, delc=12.5):

    with open(f'./TW_3a_cond_hydrofacies_realisations/hydrofacies_categories/{field_name}.grd', 'r') as file:
        ncol, nrow, nlay = map(int, file.readline().split())

        # Read the entire file into a flattened NumPy array
        iarray_flat = np.fromfile(file, dtype=float, count=-1, sep=' ')

        # Reshape the flattened array into a 3D array. Fortran-style ordering (column-major)
        iarray = iarray_flat.reshape((ncol, nrow, nlay), order='F')
        # reshape for flopy
        iarray = np.abs(iarray.transpose((2, 1, 0)))
        # truncate
        iarray = iarray[:int((H/100)*nlay),:,:]
        iarray = iarray[:,:int((W/3600)*nrow),:]
        iarray = iarray[:,:,:int((L/4400)*ncol)]

        # import hydrofacies properties
        i = int(field_name.strip("TSim_Out")) - 1
        facies_k = np.genfromtxt('./TW_3a_cond_hydrofacies_realisations/hydrofacies_conductivity_list.csv', skip_header=1, delimiter=",")[i, 2:]
        # change from facies to k
        iarray[iarray==1] = facies_k[0] 
        iarray[iarray==2] = facies_k[1]
        iarray[iarray==3] = facies_k[2]
        # increase_frequency of sample
        iarray = iarray.repeat(int(10/delv), axis=0)
        iarray = iarray.repeat(int(12.5/delr), axis=1)
        iarray = iarray.repeat(int(12.5/delc), axis=2)
        # add aquitard
        hk = iarray
        hka = [1 for i in range(int(H/delv))] #  # [13/50 for i in range(int(H/delv))]#
        return hk, hka