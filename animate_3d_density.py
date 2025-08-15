import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.ndimage import gaussian_filter

def load_density_data(filename, Nx, Ny, norm, sigma):
    """Load empirical density data from binary file"""
    with open(filename, 'rb') as f:
        raw = f.read()

    offset = 0
    times = []
    rhos = []
    
    while offset < len(raw):
        # Read time
        t_macro = np.frombuffer(raw, np.float64, 1, offset)[0]
        offset += 8
        
        # Read empirical density
        empirical = np.frombuffer(raw, np.float64, Nx * Ny, offset).reshape(Nx, Ny)
        offset += Nx * Ny * 8
        
        # Skip flow data (4 components per site)
        offset += 4 * Nx * Ny * 8
        
        # Skip marker
        offset += 4
        
        # Apply smoothing and normalization
        rho = gaussian_filter(empirical / norm, sigma=sigma)
        
        times.append(t_macro)
        rhos.append(rho)
    
    return times, rhos

def create_3d_animation(times, rhos, output_path, fps=10):
    """Create animated 3D visualization of density evolution"""
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create grid
    X, Y = np.meshgrid(np.arange(rhos[0].shape[1]), np.arange(rhos[0].shape[0]))
    
    # Find global min/max for consistent coloring
    vmin = min(rho.min() for rho in rhos)
    vmax = max(rho.max() for rho in rhos)
    
    # Plot initial surface
    surf = ax.plot_surface(X, Y, rhos[0], cmap=cm.viridis,
                           linewidth=0, antialiased=True, 
                           rcount=100, ccount=100,
                           vmin=vmin, vmax=vmax)
    
    # Add color bar
    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='Density')
    
    # Title with time
    title = ax.set_title(f'Empirical Density at t = {times[0]:.2f}', fontsize=14)
    ax.set_xlabel('X Position', fontsize=12)
    ax.set_ylabel('Y Position', fontsize=12)
    ax.set_zlabel('Density', fontsize=12)
    ax.view_init(elev=30, azim=45)
    
    # Animation update function
    def update(frame):
        ax.clear()
        surf = ax.plot_surface(X, Y, rhos[frame], cmap=cm.viridis,
                               linewidth=0, antialiased=True, 
                               rcount=100, ccount=100,
                               vmin=vmin, vmax=vmax)
        ax.set_title(f'Empirical Density at t = {times[frame]:.2f}', fontsize=14)
        ax.set_xlabel('X Position', fontsize=12)
        ax.set_ylabel('Y Position', fontsize=12)
        ax.set_zlabel('Density', fontsize=12)
        ax.view_init(elev=30, azim=45)
        return surf,
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=len(rhos), 
                        interval=1000//fps, blit=False)
    
    # Save animation
    anim.save(output_path, writer=PillowWriter(fps=fps), dpi=150)
    plt.close()
    print(f"3D animation saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Animated 3D Visualization of Empirical Density")
    parser.add_argument('--file', required=True, help="Path to stats.bin file")
    parser.add_argument('--size', type=int, default=100, help="Grid size (NxN)")
    parser.add_argument('--particles', type=int, required=True, help="Number of particles")
    parser.add_argument('--sigma', type=float, default=2.0, help="Smoothing parameter")
    parser.add_argument('--mean', type=int, default=10, help="Mean steps")
    parser.add_argument('--output', required=True, help="Output directory")
    parser.add_argument('--fps', type=int, default=5, help="Frames per second")
    args = parser.parse_args()

    Nx = args.size
    Ny = args.size
    norm = args.mean * args.particles

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Load data
    times, rhos = load_density_data(args.file, Nx, Ny, norm, args.sigma)
    print(f"Loaded {len(times)} time points")
    
    # Create output path
    output_path = os.path.join(args.output, "density_3d_animation.gif")
    
    # Create animation
    create_3d_animation(times, rhos, output_path, args.fps)

if __name__ == '__main__':
    main()
