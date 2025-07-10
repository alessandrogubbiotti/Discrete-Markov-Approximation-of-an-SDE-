import numpy as np
import struct
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

def load_data(filename, N_particles, T_macro):
    int_struct = struct.Struct('i')      # For time and marker
    coords_struct = struct.Struct('ii')  # For (x, y) particle positions

    records = []

    with open(filename, 'rb') as f:
        timestep = 0
        while True:
            time_data = f.read(int_struct.size)
            if len(time_data) < int_struct.size:
                break
            (time_raw,) = int_struct.unpack(time_data)

            particles = []
            for pid in range(N_particles):
                pdata = f.read(coords_struct.size)
                if len(pdata) < coords_struct.size:
                    raise EOFError("Incomplete particle data")
                x, y = coords_struct.unpack(pdata)
                particles.append((pid, x, y))

            marker_data = f.read(int_struct.size)
            if len(marker_data) < int_struct.size:
                raise EOFError("Missing marker")
            (marker,) = int_struct.unpack(marker_data)
            if marker != 60:
                raise ValueError(f"Expected marker 60, got {marker} at timestep {timestep}")

            for pid, x, y in particles:
                records.append((time_raw, pid, x, y))
            timestep += 1

    # Convert to structured numpy array
    dtype = np.dtype([('time', 'i4'), ('particle_id', 'i4'), ('x', 'i4'), ('y', 'i4')])
    data = np.array(records, dtype=dtype)

    # Normalize time
    n_timesteps = len(data) // N_particles
    data['time'] = data['time'] / n_timesteps * T_macro

    return data, N_particles, n_timesteps

def empirical_measure(data, lattice_size, n_timesteps, n_particles):
    density = np.zeros((lattice_size, lattice_size))
    for entry in data:
        x, y = entry['x'] % lattice_size, entry['y'] % lattice_size
        density[x, y] += 1
    return density / (n_timesteps * n_particles)

def empirical_flow(data, lattice_size, n_particles):
    flow_x = np.zeros((lattice_size, lattice_size))
    flow_y = np.zeros((lattice_size, lattice_size))

    particle_data = defaultdict(list)
    for d in data:
        particle_data[d['particle_id']].append((d['time'], d['x'], d['y']))

    for traj in particle_data.values():
        traj.sort(key=lambda x: x[0])
        for i in range(1, len(traj)):
            _, x0, y0 = traj[i - 1]
            _, x1, y1 = traj[i]

            dx = (x1 - x0 + lattice_size) % lattice_size
            dy = (y1 - y0 + lattice_size) % lattice_size
            if dx == lattice_size - 1:
                dx = -1
            elif dx >= 2:
                dx = 0
            if dy == lattice_size - 1:
                dy = -1
            elif dy >= 2:
                dy = 0

            flow_x[x0 % lattice_size, y0 % lattice_size] += dx
            flow_y[x0 % lattice_size, y0 % lattice_size] += dy

    return flow_x / (n_particles * (lattice_size ** 2)), flow_y / (n_particles * (lattice_size ** 2))

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
    plt.ylim(0, flow_y.shape[1])
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.show()

def plot_current_positions(data, lattice_size):
    last_time = data['time'].max()
    final_positions = data[data['time'] == last_time]

    plt.figure(figsize=(6, 6))
    plt.scatter(final_positions['x'] % lattice_size, final_positions['y'] % lattice_size, c='blue', s=10)
    plt.title(f"Particle Positions at Final Time (t={last_time:.2f})")
    plt.xlim(0, lattice_size)
    plt.ylim(0, lattice_size)
    plt.gca().set_aspect('equal')
    plt.grid(True)
    plt.show()

def plot_trajectories(data, lattice_size):
    plt.figure(figsize=(6, 6))
    for pid in np.unique(data['particle_id']):
        traj = data[data['particle_id'] == pid]
        plt.plot(traj['x'] % lattice_size, traj['y'] % lattice_size, alpha=0.3)
    plt.title("Particle Trajectories")
    plt.xlim(0, lattice_size)
    plt.ylim(0, lattice_size)
    plt.gca().set_aspect('equal')
    plt.grid(True)
    plt.show()
import argparse
import numpy as np

# (Assumes all previously defined functions from earlier cell are included here: 
# load_data, empirical_measure, empirical_flow, plot_empirical_density, plot_empirical_flow, 
# plot_current_positions, plot_trajectories)

def main():
    parser = argparse.ArgumentParser(description="Analyze particle trajectories from binary simulation.")
    parser.add_argument("file", type=str, help="Path to the binary trajectory file (e.g., trajectory.bin)")
    parser.add_argument("--N_particles", type=int, required=True, help="Number of particles simulated")
    parser.add_argument("--lattice_size", type=int, required=True, help="Size of the lattice (NxN)")
    parser.add_argument("--T", type=float, default=1.0, help="Macroscopic final time (default: 1.0)")

    args = parser.parse_args()

    print("Loading data...")
    data, n_particles, n_timesteps = load_data(args.file, args.N_particles, args.T)
    print(f"Loaded {len(data)} records for {n_particles} particles over {n_timesteps} timesteps.")

    print("Computing empirical measure...")
    density = empirical_measure(data, args.lattice_size, n_timesteps, n_particles)

    print("Computing empirical flow...")
    flow_x, flow_y = empirical_flow(data, args.lattice_size, n_particles)

    print("Plotting empirical density...")
    plot_empirical_density(density)

    print("Plotting empirical flow...")
    plot_empirical_flow(flow_x, flow_y)

    print("Plotting final particle positions...")
    plot_current_positions(data, args.lattice_size)

    print("Plotting particle trajectories...")
    plot_trajectories(data, args.lattice_size)

if __name__ == "__main__":
    main()

