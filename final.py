import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import struct

def load_data(filename, N_particles, T_macro):
    int_size = 4  # 4 bytes per int
    record_size = 1 + 2 * N_particles + 1  # time + particles + marker
    dtype = np.dtype([('time', 'i4'),
                      ('x', f'({N_particles},)i4'),
                      ('y', f'({N_particles},)i4')])

    raw = np.fromfile(filename, dtype=np.int32)
    
    if len(raw) % record_size != 0:
        raise ValueError("Invalid file size — does not match expected structure.")

    n_timesteps = len(raw) // record_size
    data = raw.reshape((n_timesteps, record_size))

    # Check markers
    markers = data[:, -1]
    if not np.all(markers == 60):
        idx = np.where(markers != 60)[0][0]
        raise ValueError(f"Expected marker 60, got {markers[idx]} at timestep {idx}")

    times = data[:, 0]
    coords = data[:, 1:-1].reshape(n_timesteps, N_particles, 2)
    x = coords[:, :, 0]
    y = coords[:, :, 1]

    # Normalize time to macroscopic time
    times = times / n_timesteps * T_macro

    return times, x, y, n_timesteps

def empirical_measure(x, y, N, n_timesteps, n_particles):
    x_flat = x % N
    y_flat = y % N
    density = np.zeros((N, N), dtype=np.float64)
    
    for i in range(n_timesteps):
        np.add.at(density, (x_flat[i], y_flat[i]), 1)

    return density / (n_timesteps * n_particles)

def empirical_flow(x, y, N):
    dx = (x[1:] - x[:-1]) % N
    dy = (y[1:] - y[:-1]) % N

    # Handle periodic boundaries
    dx[dx == N - 1] = -1
    dx[dx >= 2] = 0
    dy[dy == N - 1] = -1
    dy[dy >= 2] = 0

    x0 = x[:-1] % N
    y0 = y[:-1] % N

    flow_x = np.zeros((N, N))
    flow_y = np.zeros((N, N))

    for i in range(dx.shape[0]):
        np.add.at(flow_x, (x0[i], y0[i]), dx[i])
        np.add.at(flow_y, (x0[i], y0[i]), dy[i])

    norm = x.shape[0] * x.shape[1]  # total number of particle steps
    return flow_x / norm, flow_y / norm

import matplotlib.animation as animation
import matplotlib.animation as animation

def animate_particles(x, y, times, N, interval=100):
    """
    Show an animated scatter plot of particles.
    Frames are shown once every N^2 timesteps.
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    scat = ax.scatter([], [], s=10, c='blue')
    ax.set_xlim(0, N)
    ax.set_ylim(0, N)
    ax.set_title("Particle Animation")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect('equal')
    ax.grid(True)

    step_interval = N * N
    frames = np.arange(0, len(x), step_interval)

    def init():
        scat.set_offsets([])
        return scat,

    def update(frame_idx):
        xi = x[frame_idx]
        yi = y[frame_idx]
        coords = np.stack([xi, yi], axis=1)
        scat.set_offsets(coords)
        ax.set_title(f"Particle Animation (t = {times[frame_idx]:.2f})")
        return scat,

    ani = animation.FuncAnimation(
        fig, update, frames=frames,
        init_func=init, blit=True, interval=interval
    )

    plt.tight_layout()
    plt.show()


def plot_empirical_density(density):
    plt.figure(figsize=(6, 6))
    sns.heatmap(density.T, cmap='viridis', square=True, cbar_kws={"label": "Empirical Density"})
    plt.title("Empirical Measure (Heatmap)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.show()

def plot_empirical_flow(flow_x, flow_y):
    plt.figure(figsize=(6, 6))
    X, Y = np.meshgrid(np.arange(flow_x.shape[0]), np.arange(flow_x.shape[1]), indexing='ij')
    plt.quiver(X, Y, flow_x, flow_y, color='red')
    plt.title("Empirical Flow (Quiver)")
    plt.xlim(0, flow_x.shape[0])
    plt.ylim(0, flow_x.shape[1])
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.show()

def plot_current_positions(x, y, times, N):
    final_x = x[-1]
    final_y = y[-1]

    plt.figure(figsize=(6, 6))
    plt.scatter(final_x, final_y, c='blue', s=10)
    plt.title(f"Particle Positions at Final Time (t={times[-1]:.2f})")
    plt.xlim(0, N)
    plt.ylim(0, N)
    plt.gca().set_aspect('equal')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process and plot particle data.")
    parser.add_argument("--file", type=str, required=True, help="Path to binary trajectory file")
    parser.add_argument("--particles", type=int, required=True, help="Number of particles")
    parser.add_argument("--size", type=int, required=True, help="Lattice size (NxN)")
    parser.add_argument("--T", type=float, default=1.0, help="Final macroscopic time")
    args = parser.parse_args()

    print("Loading data...")
    times, x, y, n_timesteps = load_data(args.file, args.particles, args.T)
    print(f"Loaded data: {n_timesteps} timesteps, {args.particles} particles")

    print("Computing empirical measure...")
    density = empirical_measure(x, y, args.size, n_timesteps, args.particles)

    print("Computing empirical flow...")
    flow_x, flow_y = empirical_flow(x, y, args.size)

    print("Displaying particle animation...")
    animate_particles(x, y, times, args.size)

    plot_empirical_density(density)
    plot_empirical_flow(flow_x, flow_y)
    plot_current_positions(x, y, times, args.size)



