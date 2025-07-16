import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers
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
        t = int.from_bytes(raw[offset:offset + 4], 'little')
        offset += 4

        empirical = np.frombuffer(raw, np.float64, Nx * Ny, offset).reshape(Nx, Ny)
        offset += Nx * Ny * 8
        north = np.frombuffer(raw, np.float64, Nx * Ny, offset).reshape(Nx, Ny)
        offset += Nx * Ny * 8
        south = np.frombuffer(raw, np.float64, Nx * Ny, offset).reshape(Nx, Ny)
        offset += Nx * Ny * 8
        east = np.frombuffer(raw, np.float64, Nx * Ny, offset).reshape(Nx, Ny)
        offset += Nx * Ny * 8
        west = np.frombuffer(raw, np.float64, Nx * Ny, offset).reshape(Nx, Ny)
        offset += Nx * Ny * 8
        offset += 4  # skip marker

        rho = gaussian_filter(empirical / norm, sigma=sigma)
        north_s = gaussian_filter(north, sigma=sigma)
        south_s = gaussian_filter(south, sigma=sigma)
        east_s = gaussian_filter(east, sigma=sigma)
        west_s = gaussian_filter(west, sigma=sigma)

        jx = (east_s - west_s) * time_step
        jy = (north_s - south_s) * time_step

        times.append(t)
        rhos.append(rho)
        jxs.append(jx)
        jys.append(jy)

    return times, rhos, jxs, jys


def animate_field(times, rhos, jxs, jys, skip, output_folder=None):
    fig, ax = plt.subplots()
    X, Y = np.meshgrid(np.arange(rhos[0].shape[1]), np.arange(rhos[0].shape[0]))
    vmax = max(rho.max() for rho in rhos)

    im = ax.imshow(rhos[0], origin='lower', cmap='hot', vmin=0, vmax=vmax)
    Q = ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                  jxs[0][::skip, ::skip], jys[0][::skip, ::skip],
                  color='cyan', scale_units='xy', angles='xy', scale=0.005)
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, color='white', fontsize=12)

    def update(frame):
        im.set_data(rhos[frame])
        Q.set_UVC(jxs[frame][::skip, ::skip], jys[frame][::skip, ::skip])
        time_text.set_text(f't = {times[frame]}')
        return im, Q, time_text

    ax.set_title('Empirical Density and Current')
    plt.colorbar(im, ax=ax, label='Empirical density')

    anim = FuncAnimation(fig, update, frames=len(rhos), interval=200)

    # Show the animation live
    plt.show()


#def animate_field(times, rhos, jxs, jys, skip, output_folder):
#    fig, ax = plt.subplots()
#    X, Y = np.meshgrid(np.arange(rhos[0].shape[1]), np.arange(rhos[0].shape[0]))
#    vmax = max(rho.max() for rho in rhos)
#
#    im = ax.imshow(rhos[0], origin='lower', cmap='hot', vmin=0, vmax=vmax)
#    Q = ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
#                  jxs[0][::skip, ::skip], jys[0][::skip, ::skip],
#                  color='cyan', scale_units='xy', angles='xy', scale=0.005)
#    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, color='white', fontsize=12)
#
#    def update(frame):
#        im.set_data(rhos[frame])
#        Q.set_UVC(jxs[frame][::skip, ::skip], jys[frame][::skip, ::skip])
#        time_text.set_text(f't = {times[frame]}')
#        return im, Q, time_text
#
#    ax.set_title('Empirical Density and Current')
#    plt.colorbar(im, ax=ax, label='Empirical density')
#    anim = FuncAnimation(fig, update, frames=len(rhos), interval=200)
#
#    video_path = os.path.join(output_folder, 'density_and_current.gif')
#    anim.save(video_path, writer='pillow')
#    print(f"Animation saved to {video_path}")
#    plt.close()
#

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
    parser.add_argument('--debug_time', type=int, default=None)
    args = parser.parse_args()

    Nx = args.size
    Ny = args.size
    norm = args.mean * args.particles
    time_step = args.Tmacro / (args.mean * args.particles)

    output_folder = f"SDE_Potential={args.potential_type}_size={Nx}_N_particles={args.particles}"
    os.makedirs(output_folder, exist_ok=True)

    times, rhos, jxs, jys = load_data(args.file, Nx, Ny, norm, time_step, args.sigma)

    if args.debug_time is not None and args.debug_time in times:
        idx = times.index(args.debug_time)
        rho = rhos[idx]
        jx = jxs[idx]
        jy = jys[idx]
        print(f"DEBUG: t={args.debug_time}")
        print(f"Empirical measure: min={rho.min():.5f}, max={rho.max():.5f}, mean={rho.mean():.5f}")
        print(f"Current jx: min={jx.min():.5f}, max={jx.max():.5f}, mean={jx.mean():.5f}")
        print(f"Current jy: min={jy.min():.5f}, max={jy.max():.5f}, mean={jy.mean():.5f}")

        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(rho, origin='lower', cmap='hot')
        ax.quiver(np.arange(0, Nx, args.skip), np.arange(0, Ny, args.skip),
                  jx[::args.skip, ::args.skip], jy[::args.skip, ::args.skip],
                  color='cyan', scale_units='xy', angles='xy', scale=0.005)
        ax.set_title(f"Debug Plot at t={args.debug_time}")
        plt.colorbar(im, ax=ax, label='Empirical density')
        debug_plot_path = os.path.join(output_folder, f"debug_plot_t{args.debug_time}.png")
        plt.savefig(debug_plot_path)
        print(f"Debug plot saved to {debug_plot_path}")
        plt.close()

    animate_field(times, rhos, jxs, jys, args.skip, output_folder)


if __name__ == '__main__':
    main()
