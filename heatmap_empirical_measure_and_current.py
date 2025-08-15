import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.ndimage import gaussian_filter

def load_data(filename, Nx, Ny, norm, time_step, sigma):
    with open(filename, 'rb') as f:
        raw = f.read()

    offset = 0
    times = []
    rhos = []
    jxs = []
    jys = []
    
    while offset < len(raw):
        t_macro = np.frombuffer(raw, np.float64, 1, offset)[0]
        offset += 8

        empirical = np.frombuffer(raw, np.float64, Nx * Ny, offset).reshape(Nx, Ny)
        offset += Nx * Ny * 8

        # Read flow data as array of structs and reshape
        flow_data = np.frombuffer(raw, np.float64, 4 * Nx * Ny, offset).reshape(Nx, Ny, 4)
        offset += 4 * Nx * Ny * 8
        
        # Extract components from the struct
        north = flow_data[:, :, 0]
        south = flow_data[:, :, 1]
        east = flow_data[:, :, 2]
        west = flow_data[:, :, 3]

        offset += 4  # skip marker

        # Apply smoothing
        rho = gaussian_filter(empirical / norm, sigma=sigma)
        north_s = gaussian_filter(north, sigma=sigma)
        south_s = gaussian_filter(south, sigma=sigma)
        east_s = gaussian_filter(east, sigma=sigma)
        west_s = gaussian_filter(west, sigma=sigma)

        # Compute currents
        jx = 10 * (east_s - west_s) * time_step
        jy = 10 * (north_s - south_s) * time_step

        times.append(t_macro)
        rhos.append(rho)
        jxs.append(jx)
        jys.append(jy)

    return times, rhos, jxs, jys

def animate_field(times, rhos, jxs, jys, skip, output_folder=None):
    fig, ax = plt.subplots()
    vmax = max(rho.max() for rho in rhos)

    # Create coordinate grid with proper indexing
    X, Y = np.meshgrid(np.arange(rhos[0].shape[1]), 
                       np.arange(rhos[0].shape[0]),
                       indexing='xy')
    
    im = ax.imshow(rhos[0], origin='upper', cmap='hot', vmin=0, vmax=vmax)
    
    # Quiver plot with proper coordinates
    Q = ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                  jxs[0][::skip, ::skip], jys[0][::skip, ::skip],
                  color='cyan', scale_units='xy', angles='xy', scale=0.005)
    
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, color='white', fontsize=12)

    def update(frame):
        im.set_data(rhos[frame])
        Q.set_UVC(jxs[frame][::skip, ::skip], jys[frame][::skip, ::skip])
        time_text.set_text(f't = {times[frame]:.2f}')
        return im, Q, time_text

    ax.set_title('Empirical Density and Current')
    plt.colorbar(im, ax=ax, label='Empirical density')
    anim = FuncAnimation(fig, update, frames=len(rhos), interval=200)
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Visualize empirical density & current")
    parser.add_argument('--file', required=True)
    parser.add_argument('--particles', type=int, required=True)
    parser.add_argument('--size', type=int, default=100)
    parser.add_argument('--sigma', type=float, default=0.05)
    parser.add_argument('--Tmacro', type=float, default=1.0)
    parser.add_argument('--mean', type=int, default=10)
    parser.add_argument('--skip', type=int, default=5)
    parser.add_argument('--potential_type', type=str, required=True)
    parser.add_argument('--debug_time', type=float, default=None)
    args = parser.parse_args()

    Nx = args.size
    Ny = args.size
    norm = args.mean * args.particles
    time_step = args.Tmacro / (args.mean * args.particles)

    output_folder = f"SDE_Potential={args.potential_type}_size={Nx}_N_particles={args.particles}"
    os.makedirs(output_folder, exist_ok=True)

    times, rhos, jxs, jys = load_data(args.file, Nx, Ny, norm, time_step, args.sigma)

    if args.debug_time is not None:
        # Find closest time index
        idx = np.argmin(np.abs(np.array(times) - args.debug_time))
        debug_time = times[idx]
        
        print(f"DEBUG: t={debug_time:.2f} (closest to {args.debug_time})")
        print(f"Empirical measure: min={rhos[idx].min():.5f}, max={rhos[idx].max():.5f}, mean={rhos[idx].mean():.5f}")
        print(f"Current jx: min={jxs[idx].min():.5f}, max={jxs[idx].max():.5f}, mean={jxs[idx].mean():.5f}")
        print(f"Current jy: min={jys[idx].min():.5f}, max={jys[idx].max():.5f}, mean={jys[idx].mean():.5f}")

        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(rhos[idx], origin='upper', cmap='hot')
        
        # Proper coordinate handling for debug plot
        x = np.arange(0, Ny)
        y = np.arange(0, Nx)
        X, Y = np.meshgrid(x, y, indexing='xy')
        
        ax.quiver(X[::args.skip, ::args.skip], 
                  Y[::args.skip, ::args.skip],
                  jxs[idx][::args.skip, ::args.skip], 
                  jys[idx][::args.skip, ::args.skip],
                  color='cyan', scale_units='xy', angles='xy', scale=0.005)
        
        ax.set_title(f"Debug Plot at t={debug_time:.2f}")
        plt.colorbar(im, ax=ax, label='Empirical density')
        debug_plot_path = os.path.join(output_folder, f"debug_plot_t{debug_time:.2f}.png".replace('.', '_'))
        plt.savefig(debug_plot_path)
        print(f"Debug plot saved to {debug_plot_path}")
        plt.close()

    animate_field(times, rhos, jxs, jys, args.skip, output_folder)

if __name__ == '__main__':
    main()


