import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable

alpha = 0.14 # Steak diffusivity (mm^2/s)
L = 30 # Steak thickness (mm)
Ny = 50
dy = L/Ny
dt = 0.9 * dy**2 / (2*alpha)
A = alpha * dt / dy**2

T_raw = 23.0
T_air = 23.0
T_myosin = 40.0
T_glycogen = 55.0
T_myoglobin = 60.0
T_actin = 70.0
T_browning = 120.0
T_charring = 180.0


# recipe format is [time, temp_bot, temp_top]
recipe = [[10, 230, T_air], # sear bot
          [10, T_air, 230], # sear top
          [360, 110, 110], # bake
          [325, T_air, T_air]] # rest

TB = np.hstack([np.ones(int(step[0]/dt)+1)*step[1] for step in recipe])
TT = np.hstack([np.ones(int(step[0]/dt)+1)*step[2] for step in recipe])
Nt = len(TB)

# Temperature profile
T = np.ones((Nt,Ny)) * T_raw
T[:, 0] = TB
T[:, -1] = TT

# Protein profile
P = np.ones((Nt, Ny)) * T_raw

# Diffuse
for n in range(Nt-1):
    T[n+1, 1:-1] = T[n, 1:-1] + A * (T[n, 2:] - 2*T[n, 1:-1] + T[n, :-2])
    P[n+1] = np.maximum(P[n], T[n+1])

## Create Plots
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, constrained_layout=True)

t = np.arange(Nt) * dt
y = np.linspace(0, L, Ny)

# Remove the heat sources
y_plot = y[1:-1]
T_plot = T[:, 1:-1]
P_plot = P[:, 1:-1]

# Temperature heatmap
im = ax1.pcolormesh(t, y_plot, T_plot.T, cmap="turbo",
                    norm = mc.Normalize(vmin=0, vmax=100)
                    )
ax1.set_title("Heat Map of Steak")
ax1.set_ylabel("Depth (mm)")
fig.colorbar(im, ax=ax1, pad=-.2, label="Temperature (°C)")


# Define steak colormap
p_colors = {"Raw":"darkred", "Rare":"red", "Med Rare":"tomato", "Medium":"salmon", "Well":"peru", "Carmelized":"saddlebrown", "Burnt":"black"}
patches = [mpatches.Patch(color=col, label=lab) for lab, col in p_colors.items()]
steak_cmap = mc.ListedColormap(p_colors.values())
steak_norm = mc.BoundaryNorm([0, T_myosin, T_glycogen, T_myoglobin, T_actin, T_browning, T_charring,5000], ncolors=steak_cmap.N)

# Protein state map
im2 = ax2.pcolormesh(t, y_plot, P_plot.T, cmap=steak_cmap, norm=steak_norm)

ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Depth (mm)')
ax2.set_title('Protein State Map')
ax2.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
plt.savefig("steak.png", dpi=300, bbox_inches="tight")
plt.show()



































































































