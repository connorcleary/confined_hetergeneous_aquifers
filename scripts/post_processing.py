import numpy as np

def hf(conc, head, nlay, z_b, dz):
    # calculate hf using method outlined by simmons and post in the book variable density flow
    hf = np.zeros_like(conc)
    for col in range(conc.shape[1]):
        for row in range(conc.shape[0]):
            row = nlay-row-1
            zi = (z_b-dz/2-dz*row) 
            # print(row)
            hf[row, col] = zi + (head[row, col]-zi)*(1000+0.7143*conc[row,col])/1000
            # print((1000+0.7143*conc[row,col])/1000)
    return hf

def find_max_toe_position(conc, delc, x_onshore):
    # find the toe position
    interface = [-x_onshore + delc*np.min(np.argmax([conc[:, row, :] >= 0.35], axis=2)) for row in range(conc.shape[1])]
    max = np.min(interface)
    return max

def find_average_toe_position(conc, delc, x_onshore):
    interface = [-x_onshore + delc*np.min(np.argmax([conc[:, row, :] >= 0.35], axis=2)) for row in range(conc.shape[1])]
    ave = np.mean(interface)
    return ave

def find_average_width(conc, delc, ncol, x_onshore):
    extent = np.argmax(conc >= 31.5, axis=2)
    extent[extent == 0] = ncol
    return delc * np.mean(extent - np.argmax(conc >= 3.5, axis=2))

def find_toe_position_2D():
    pass

def find_width_2D():
    pass
