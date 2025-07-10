import argparse


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
import numpy as np
import struct

def load_trajectory_data(filename, N_particles, N, T_macro, stride):
    int_struct = struct.Struct('i')
    coords_struct = struct.Struct('ii')

    record_size = 4 + N_particles * 8 + 4  # time + coords + marker

    x_all = []
    y_all = []
    times = []

    with open(filename, 'rb') as f:
        frame = 0
        while True:
            time_data = f.read(int_struct.size)
            if len(time_data) < int_struct.size:
                break
            (t_raw,) = int_struct.unpack(time_data)

            coords_data = f.read(coords_struct.size * N_particles)
            if len(coords_data) < coords_struct.size * N_particles:
                break

            coords = list(coords_struct.iter_unpack(coords_data))
            xs, ys = zip(*coords)

            marker_data = f.read(int_struct.size)
            if len(marker_data) < int_struct.size:
                break
            (marker,) = int_struct.unpack(marker_data)
            if marker != 60:
                raise ValueError(f"Expected marker 60, got {marker} at frame {frame}")

            if frame % stride == 0:
                x_all.append(xs)
                y_all.append(ys)
                times.append(t_raw)

            frame += 1

    x_arr = np.array(x_all, dtype=np.int32)
    y_arr = np.array(y_all, dtype=np.int32)
    times = np.array(times, dtype=np.float32)
    times = times / times[-1] * T_macro  # Normalize to macro time

    return x_arr, y_arr, times

def compute_density_hist(x_frame, y_frame, N, N_particles):
    """Compute normalized density using vectorized histogram."""
    H, _, _ = np.histogram2d(x_frame, y_frame, bins=N, range=[[0, N], [0, N]])
    return H.T / N_particles  # transpose to match (x, y) axes for imshow

def plot_height_function(x_frame, y_frame, N, N_particles):
    """Plot the empirical (normalized) height function at a single frame."""
    height = compute_density_hist(x_frame, y_frame, N, N_particles)
    plt.figure(figsize=(6, 6))
    sns.heatmap(height, cmap="viridis", square=True, cbar_kws={"label": "Normalized Height"})
    plt.title("Height Function (Particle Density)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.show()


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns

def animate_density(x, y, times, N, n_particles, stride=200):
    """
    Animate the evolving density heatmap from particle positions.
    
    Args:
      x, y: arrays of shape (frames, n_particles) with positions
      times: array of shape (frames,) with times
      N: lattice size (NxN)
      n_particles: total number of particles
      stride: steps between frames to animate
    
    """
    # Precompute densities for selected frames only
    frames_idx = np.arange(0, len(times), stride)
    densities = []
    for idx in frames_idx:
        density = np.zeros((N, N))
        xi, yi = x[idx], y[idx]
        np.add.at(density, (xi, yi), 1)
        densities.append(density / n_particles)
    densities = np.array(densities)
    
    fig, ax = plt.subplots(figsize=(6,6))
    sns.set_style("white")
    heatmap = ax.imshow(densities[0], cmap='viridis', vmin=0, vmax=densities.max())
    cbar = plt.colorbar(heatmap, ax=ax)
    cbar.set_label('Normalized Particle Density')
    title = ax.set_title(f"Density at t = {times[frames_idx[0]]:.2f}")
    
    def update(frame):
        heatmap.set_data(densities[frame])
        title.set_text(f"Density at t = {times[frames_idx[frame]]:.2f}")
        return heatmap, title
    
    ani = animation.FuncAnimation(
        fig, update, frames=len(frames_idx),
        interval=10, blit=True, repeat=False
    )
    
    plt.show()

# Usage example (assuming you have loaded x, y, times, N, n_particles):
# animate_density(x, y, times, N=100, n_particles=1000, stride=200)

#def animate_density(x, y, times, N, N_particles, interval=200):
#    """Animate the normalized density heatmap over time."""
#    fig, ax = plt.subplots(figsize=(6, 6))
#    density = compute_density_hist(x[0], y[0], N, N_particles)
#    im = ax.imshow(density, cmap='viridis', origin='lower', vmin=0, vmax=1)
#    cb = plt.colorbar(im, ax=ax, label="Normalized Density")
#    ax.set_title(f"Density at t = {times[0]:.2f}")
#    ax.set_xlabel("x")
#    ax.set_ylabel("y")
#
#    def update(i):
#        density = compute_density_hist(x[i], y[i], N, N_particles)
#        im.set_data(density)
#        ax.set_title(f"Density at t = {times[i]:.2f}")
#        return [im]
#
#    ani = animation.FuncAnimation(fig, update, frames=len(times), interval=interval, blit=True)
#    plt.tight_layout()
#    plt.show()

def animate_log_density(x, y, times, N, N_particles, interval=200):
    """Animate the log of the normalized density heatmap."""
    fig, ax = plt.subplots(figsize=(6, 6))
    density = compute_density_hist(x[0], y[0], N, N_particles)
    log_density = np.log(density + 1e-10)  # small epsilon to avoid log(0)
    im = ax.imshow(log_density, cmap='plasma', origin='lower', vmin=np.min(log_density), vmax=0)
    cb = plt.colorbar(im, ax=ax, label="log(Density)")
    ax.set_title(f"log(Density) at t = {times[0]:.2f}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    def update(i):
        density = compute_density_hist(x[i], y[i], N, N_particles)
        log_density = np.log(density + 1e-10)
        im.set_data(log_density)
        ax.set_title(f"log(Density) at t = {times[i]:.2f}")
        return [im]

    ani = animation.FuncAnimation(fig, update, frames=len(times), interval= 10 , blit=True)
    plt.tight_layout()
    plt.show()



def main():
    parser = argparse.ArgumentParser(description="Visualize particle density evolution.")
    parser.add_argument("--file", type=str, required=True, help="Path to binary trajectory file")
    parser.add_argument("--particles", type=int, required=True, help="Number of particles")
    parser.add_argument("--size", type=int, required=True, help="Grid size (NxN)")
    parser.add_argument("--T", type=float, default=1.0, help="Macroscopic final time")
    parser.add_argument("--stride", type=int, default=200, help="Frame stride (simulation steps per saved frame)")
    args = parser.parse_args()

    print("Loading data...")
    x, y, times = load_trajectory_data(args.file, args.particles, args.size, args.T, args.stride)
    print(f"Loaded data for {x.shape[0]} timesteps with {x.shape[1]} particles each.")

    print("Plotting final height function...")
    plot_height_function(x[-1], y[-1], args.size, args.particles)

    print("Animating density heatmap...")
    animate_density(x, y, times, args.size, args.particles)

    print("Animating log-density heatmap...")
    animate_log_density(x, y, times, args.size, args.particles)

if __name__ == "__main__":
    main()

