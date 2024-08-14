import matplotlib.pyplot as plt
import numpy as np
from random_fields import import_field
import matplotlib.cm as cm
from matplotlib.colors import LogNorm, LightSource
import matplotlib as mpl
import matplotlib.patches as mpatches

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

H = 20
D = 20 
kb = 0.1
field = import_field("TSim_Out5", H, 10, D, kb, W=1000)
x_onshore = 400
x_offshore = 4000
z_slice = 30
x_slice = 2000
width = 1000
y_slice = 500

ls = LightSource(270, 35)
field = field[0][:int((H+D)/10), :int(width/12.5) ,:]
norm = LogNorm(vmin=0.1, vmax=np.max(field))
#plot left base
x = np.linspace(-x_onshore, x_offshore, int((x_onshore+x_offshore)/12.5)+1)
z = np.linspace(-z_slice, -H-D, int((H+D-z_slice)/10)+1)
X, Z = np.meshgrid(x, z)
k = field[int(z_slice/10):, 0, :]
colors = cm.coolwarm(norm(k))
ax.plot_surface(X, 0, Z, facecolors=colors, lightsource=ls)

# plot left upper 
x = np.linspace(-x_onshore, x_slice, int((x_onshore+x_slice)/12.5)+1)
z = np.linspace(0, -z_slice, int(z_slice/10)+1)
X, Z = np.meshgrid(x, z)
k = field[:int(z_slice/10), 0, :int((x_onshore+x_slice)/12.5)]
colors = cm.coolwarm(norm(k))
ax.plot_surface(X, 0, Z, facecolors=colors, lightsource=ls)

# plot top left
x = np.linspace(-x_onshore, x_slice, int((x_onshore+x_slice)/12.5)+1)
y = np.linspace(0, width, int(width/12.5)+1)
X, Y = np.meshgrid(x, y)
k = field[0, :, :int((x_onshore+x_slice)/12.5)]
colors = cm.coolwarm(norm(k))
Z = np.zeros_like(Y)
ax.plot_surface(X, Y, Z, facecolors=colors, lightsource=ls)

# plot top right
x = np.linspace(x_slice, x_offshore, int((x_offshore-x_slice)/12.5)+1)
y = np.linspace(y_slice, width, int((width-y_slice)/12.5)+1)
X, Y = np.meshgrid(x, y)
k = field[0, int(y_slice/12.5):, int((x_onshore+x_slice)/12.5):]
colors = cm.coolwarm(norm(k))
Z = np.zeros_like(Y)
ax.plot_surface(X, Y, Z, facecolors=colors, lightsource=ls)

# plot middle left
y = np.linspace(0, y_slice, int(y_slice/12.5)+1)
z = np.linspace(0, -z_slice, int(z_slice/10)+1)
Y, Z = np.meshgrid(y, z)
k = field[:int(z_slice/10), :int(y_slice/12.5), int((x_onshore+x_slice)/12.5)-2]
colors = cm.coolwarm(norm(k))
ax.plot_surface(x_slice, Y, Z, facecolors=colors, lightsource=ls)

# plot middle right
x = np.linspace(x_slice, x_offshore, int((x_offshore-x_slice)/12.5)+1)
z = np.linspace(0, -z_slice, int(z_slice/10)+1)
X, Z = np.meshgrid(x, z)
k = field[:int(z_slice/10), int(y_slice/12.5), int((x_onshore+x_slice)/12.5):]
colors = cm.coolwarm(norm(k))
ax.plot_surface(X, y_slice, Z, facecolors=colors, lightsource=ls)

# plot middle bottom
x = np.linspace(x_slice, x_offshore, int((x_offshore-x_slice)/12.5)+1)
y = np.linspace(0, y_slice, int(y_slice/12.5)+1)
X, Y = np.meshgrid(x, y)
k = field[int(z_slice/12.5), :int(y_slice/12.5), int((x_onshore+x_slice)/12.5):]
colors = cm.coolwarm(norm(k))
Z = -(z_slice)*np.ones_like(Y)
ax.plot_surface(X, Y, Z, facecolors=colors, lightsource=ls)

# plot right base
y = np.linspace(0, width, int(width/12.5)+1)
z = np.linspace(-z_slice, -H-D, int((H+D-z_slice)/10)+1)
Y, Z = np.meshgrid(y, z)
k = field[int(z_slice/10):, :, -1]
colors = cm.coolwarm(norm(k))
ax.plot_surface(x_offshore, Y, Z, facecolors=colors, lightsource=ls)

# plot right upper
y = np.linspace(y_slice, width, int((width-y_slice)/12.5)+1)
z = np.linspace(0, -z_slice, int(z_slice/10)+1)
Y, Z = np.meshgrid(y, z)
k = field[:int(z_slice/10), int(y_slice/12.5):, -1]
colors = cm.coolwarm(norm(k))
ax.plot_surface(x_offshore, Y, Z, facecolors=colors, lightsource=ls)

ax.set_xlabel("Distance offshore [m]", labelpad=10.0)
ax.set_ylabel("Distance alongshore [m]", labelpad=5.0)
ax.set_zlabel("Depth [m]")
ax.set_box_aspect((1, 0.8, 0.2))
ax.set_zticks([0, -20, -40])
ax.set_yticks([0, 500, 1000])
ax.set_yticks([0, 1000, 2000, 3000, 4000])

# m = cm.ScalarMappable(cmap=cm.coolwarm, norm=norm)
# .set_array(field)
# plt.colorbar(m, ax=ax, location="left", label=r"C [PSU]", shrink=0.75)

facies_values = np.sort(np.unique(field), axis=None)

facies0 = mpatches.Patch(color=cm.coolwarm(norm(facies_values[3])), label='')
facies1 = mpatches.Patch(color=cm.coolwarm(norm(facies_values[2])), label='The red data')
facies2 = mpatches.Patch(color=cm.coolwarm(norm(facies_values[1])), label='The red data')
aquitard = mpatches.Patch(color=cm.coolwarm(norm(facies_values[0])), label='The red data')

blue_patch = mpatches.Patch(color='blue', label='The blue data')

plt.show()