import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 5:
    print("Usage: python plot_forms.py <binary_file> <N> <skip> <name>")
    sys.exit(1)

binary_file = sys.argv[1]
N = int(sys.argv[2])
skip = int(sys.argv[3])
label = sys.argv[4]

# Read the binary file
forms = np.fromfile(binary_file, dtype=[('x','f16'), ('y','f16')]).reshape(N, N)

# Plot quiver
X, Y = np.meshgrid(np.arange(N), np.arange(N), indexing='xy')
U = forms['x']
V = forms['y']

plt.figure(figsize=(6,6))
plt.quiver(X[::skip, ::skip], Y[::skip, ::skip],
           U[::skip, ::skip], V[::skip, ::skip])
plt.title(label)
plt.show()



