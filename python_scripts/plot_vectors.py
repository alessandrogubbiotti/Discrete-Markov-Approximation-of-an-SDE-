import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

if len(sys.argv) != 6:
    print("Usage: python plot_vectors.py <binary_file> <N> <sigma> <skip> <name>")
    sys.exit(1)


filename = sys.argv[1]
N = int(sys.argv[2])
sigma = float(sys.argv[3])
skip = int(sys.argv[4])
name = sys.argv[5]

# Read binary file
data = np.fromfile(filename, dtype=np.float64).reshape(N, N, 4)
north, south, east, west = data[:, :, 0], data[:, :, 1], data[:, :, 2], data[:, :, 3]

# Gaussian smoothing
north_s = gaussian_filter(north, sigma=sigma)
south_s = gaussian_filter(south, sigma=sigma)
east_s  = gaussian_filter(east, sigma=sigma)
west_s  = gaussian_filter(west, sigma=sigma)

# Compute currents
time_step = 0.01
jx = 10 * (east_s - west_s) * time_step
jy = 10 * (north_s - south_s) * time_step

# Apply skip
X, Y = np.meshgrid(np.arange(N), np.arange(N), indexing='xy')
X_skip = X[::skip, ::skip]
Y_skip = Y[::skip, ::skip]
jx_skip = jx[::skip, ::skip]
jy_skip = jy[::skip, ::skip]

# Plot quiver
plt.figure()
plt.quiver(X_skip, Y_skip, jx_skip, jy_skip, color='cyan', scale_units='xy', angles='xy', scale=0.005)
plt.title(f"Current vectors: {name}")
plt.show()

