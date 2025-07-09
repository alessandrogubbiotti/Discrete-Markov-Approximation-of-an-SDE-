import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def load_data(filename, N, T_macro):
    dt = np.dtype([('time', 'i4'), ('particle_id', 'i4'), ('x', 'i4'), ('y', 'i4')])
    data = np.fromfile(filename, dtype=dt)

    n_particles = np.unique(data['particle_id']).size
    n_timesteps = N**2

    data = data.copy()
    data['time'] = data['time'] / n_timesteps * T_macro  # normalize to macroscopic time

    return data, n_particles, n_timesteps

def empirical_measure(data, N, n_timesteps, n_particles):
    density = np.zeros((N, N))
    for entry in data:
        x, y = entry['x'] % N, entry['y'] % N
        density[x, y] += 1
    return density / (n_timesteps * n_particles)

def empirical_flow(data, N, n_particles):
    from collections import defaultdict

    flow_x = np.zeros((N, N))
    flow_y = np.zeros((N, N))

    particle_data = defaultdict(list)
    for d in data:
        particle_data[d['particle_id']].append((d['time'], d['x'], d['y']))

    for traj in particle_data.values():
        traj.sort()
        for i in range(1, len(traj)):
            _, x0, y0 = traj[i - 1]
            _, x1, y1 = traj[i]

            dx = (x1 - x0 + N) % N
            dy = (y1 - y0 + N) % N
            if dx == N - 1:
                dx = -1
            elif dx >= 2:
                dx = 0
            if dy == N - 1:
                dy = -1
            elif dy >= 2:
                dy = 0

            flow_x[x0 % N, y0 % N] += dx
            flow_y[x0 % N, y0 % N] += dy

    return flow_x / (n_particles * (N**2)), flow_y / (n_particles * (N**2))

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

def plot_current_positions(data, N):
    last_time = data['time'].max()
    final_positions = data[data['time'] == last_time]

    plt.figure(figsize=(6, 6))
    plt.scatter(final_positions['x'], final_positions['y'], c='blue', s=10)
    plt.title(f"Particle Positions at Final Time (t={last_time:.2f})")
    plt.xlim(0, N)
    plt.ylim(0, N)
    plt.gca().set_aspect('equal')
    plt.grid(True)
    plt.show()

def plot_trajectories(data, N):
    from collections import defaultdict

    particle_data = defaultdict(list)
    for d in data:
        particle_data[d['particle_id']].append((d['x'], d['y']))

    plt.figure(figsize=(6, 6))
    for traj in particle_data.values():
        xs, ys = [], []
        for i in range(1, len(traj)):
            x0, y0 = traj[i - 1]
            x1, y1 = traj[i]

            dx = abs(x1 - x0)
            dy = abs(y1 - y0)
            if dx > 1 or dy > 1:  # jump across torus, don't draw
                if len(xs) > 1:
                    plt.plot(xs, ys, alpha=0.5)
                xs, ys = [], []
            else:
                xs.extend([x0, x1])
                ys.extend([y0, y1])
        if len(xs) > 1:
            plt.plot(xs, ys, alpha=0.5)

    plt.title("Particle Trajectories (No Wrap Lines)")
    plt.xlim(0, N)
    plt.ylim(0, N)
    plt.gca().set_aspect('equal')
    plt.grid(True)
    plt.show()

# Main entry point
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, required=True, help="Lattice size")
    parser.add_argument("--T", type=float, required=True, help="Macroscopic time (T)")
    parser.add_argument("--file", type=str, default="traj_particles.bin", help="Binary input file")
    args = parser.parse_args()

    data, n_particles, n_timesteps = load_data(args.file, args.N, args.T)
    density = empirical_measure(data, args.N, n_timesteps, n_particles)
    flow_x, flow_y = empirical_flow(data, args.N, n_particles)

    plot_empirical_density(density)
    plot_empirical_flow(flow_x, flow_y)
    plot_current_positions(data, args.N)
    plot_trajectories(data, args.N)

