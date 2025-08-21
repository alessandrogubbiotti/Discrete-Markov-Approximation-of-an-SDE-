import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.colors as mcolors

def read_trajectory_file(filename, num_particles, grid_size):
    """Read particle trajectory data from binary file"""
    particles = []
    times = []
    marker_size = 4  # Size of the marker integer in bytes
    time_size = 8    # Size of double in bytes
    particle_size = 8  # Size of two integers (2 * 4 bytes)
    
    with open(filename, 'rb') as f:
        frame_index = 0
        while True:
            # Read time
            time_data = f.read(time_size)
            if not time_data or len(time_data) < time_size:
                break
            t_macro = np.frombuffer(time_data, dtype=np.float64)[0]
            times.append(t_macro)
            
            # Read particles
            frame_particles = []
            for _ in range(num_particles):
                particle_data = f.read(particle_size)
                if not particle_data or len(particle_data) < particle_size:
                    break
                x, y = np.frombuffer(particle_data, dtype=np.int32)
                frame_particles.append((x, y))
            
            # Skip marker
            f.read(marker_size)
            
            particles.append(frame_particles)
            frame_index += 1
            
            # Print progress every 10 frames
            if frame_index % 10 == 0:
                print(f"Read frame {frame_index} at t = {t_macro:.2f}")
    
    return times, particles

def create_particle_animation(times, particles, grid_size, output_path, fps=10):
    """Create animation of particle evolution"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Create grid background
    ax.set_xlim(0, grid_size)
    ax.set_ylim(0, grid_size)
    ax.set_xticks(np.arange(0, grid_size+1, 5))
    ax.set_yticks(np.arange(0, grid_size+1, 5))
    ax.grid(True, color='gray', linestyle='-', linewidth=0.5)
    ax.set_title(f'Particle Positions at t = {times[0]:.2f}', fontsize=14)
    ax.set_xlabel('X Position', fontsize=12)
    ax.set_ylabel('Y Position', fontsize=12)
    
    # Create scatter plot with distinct colors
    num_particles = len(particles[0])
    colors = list(mcolors.TABLEAU_COLORS.values())
    if num_particles > len(colors):
        colors = colors * (num_particles // len(colors) + 1)
    
    # Create initial scatter plot
    scat = ax.scatter([], [], s=50, alpha=0.7)
    
    # FIXED: Corrected the syntax error here
    all_x = np.zeros((len(particles), num_particles))
    all_y = np.zeros((len(particles), num_particles))
    
    for frame_idx, frame_particles in enumerate(particles):
        for particle_idx, (x, y) in enumerate(frame_particles):
            all_x[frame_idx, particle_idx] = x
            all_y[frame_idx, particle_idx] = y
    
    # Animation update function
    def update(frame):
        # Update positions
        scat.set_offsets(np.c_[all_x[frame], all_y[frame]])
        
        # Update title with time
        ax.set_title(f'Particle Positions at t = {times[frame]:.2f}')
        return scat,
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=len(particles), 
                         interval=1000//fps, blit=True)
    
    # Save animation
    anim.save(output_path, writer=PillowWriter(fps=fps), dpi=100)
    plt.close()
    print(f"Particle animation saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Particle Trajectory Visualization")
    parser.add_argument('--file', required=True, help="Path to trajectory file")
    parser.add_argument('--particles', type=int, required=True, help="Number of particles")
    parser.add_argument('--size', type=int, required=True, help="Grid size (N)")
    parser.add_argument('--output', required=True, help="Output directory")
    parser.add_argument('--fps', type=int, default=10, help="Frames per second")
    parser.add_argument('--max_frames', type=int, default=1000, help="Maximum frames to process")
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    try:
        # Read trajectory data
        print(f"Reading trajectory data from {args.file}...")
        times, particles = read_trajectory_file(args.file, args.particles, args.size)
        
        # Limit frames if requested
        if args.max_frames > 0 and len(particles) > args.max_frames:
            skip = len(particles) // args.max_frames
            particles = particles[::skip]
            times = times[::skip]
            print(f"Downsampled to {len(particles)} frames")
        
        print(f"Loaded {len(particles)} frames with {args.particles} particles each")
        
        # Create output path
        output_path = os.path.join(args.output, "particle_trajectory.gif")
        
        # Create animation
        create_particle_animation(times, particles, args.size, output_path, args.fps)
    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    main()



