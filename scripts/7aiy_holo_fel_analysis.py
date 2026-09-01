import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import os

# Force matplotlib to not use a screen
matplotlib.use('Agg')

def read_xvg_column(file_path):
    """Reads data from an XVG file and converts from nm to Angstrom."""
    data = []
    if not os.path.exists(file_path):
        print(f"Error: {file_path} not found.")
        return None
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith(('@', '#')) and line.strip():
                parts = line.split()
                if len(parts) >= 2:
                    # Convert NM to Angstrom (multiplied by 10)
                    data.append(float(parts[1]) * 10)
    return np.array(data)

def plot_publication_fel(pc1_file, pc2_file, system_name):
    print(f"Generating high-quality publication FEL for {system_name}...")

    p1 = read_xvg_column(pc1_file)
    p2 = read_xvg_column(pc2_file)

    if p1 is None or p2 is None: return

    # Sync lengths
    min_len = min(len(p1), len(p2))
    x = p1[:min_len]
    y = p2[:min_len]

    # 1. Grid Resolution (Keep it around 40-50 for that specific jagged look)
    nbins = 50
    hist, xedges, yedges = np.histogram2d(x, y, bins=nbins)

    # 2. Gibbs Free Energy Calculation
    R = 0.008314  # Gas constant in kJ/mol/K
    T = 300       # Temperature in Kelvin

    prob = hist / np.sum(hist)
    # Avoid log(0) and calculate Energy
    free_energy = -R * T * np.log(prob + 1e-10)
    free_energy -= np.min(free_energy)

    # ============================================================
    # QUANTITATIVE FEL ANALYSIS (before energy cap)
    # ============================================================

    print("Calculating FEL statistics...")

    # Save original (uncapped) free energy
    analysis_energy = free_energy.copy()

    # Overall energy statistics
    mean_energy = np.mean(analysis_energy)
    max_energy = np.max(analysis_energy)
    std_energy = np.std(analysis_energy)

    # Bin centres
    x_centers = (xedges[:-1] + xedges[1:]) / 2
    y_centers = (yedges[:-1] + yedges[1:]) / 2

    # Global minimum
    min_index = np.nanargmin(analysis_energy)
    i_min, j_min = np.unravel_index(min_index, analysis_energy.shape)

    global_min_energy = analysis_energy[i_min, j_min]
    global_pc1 = x_centers[i_min]
    global_pc2 = y_centers[j_min]

    # PC ranges
    pc1_min = np.min(x)
    pc1_max = np.max(x)

    pc2_min = np.min(y)
    pc2_max = np.max(y)

    pc1_range = pc1_max - pc1_min
    pc2_range = pc2_max - pc2_min

    # 3. Setting an energy cap for cleaner visualization (like in your target image)
    # Typically, anything above 7-10 kJ/mol is shown as a flat plateau
    energy_cap = 8.0
    free_energy[free_energy > energy_cap] = energy_cap

    # Create meshgrid (reuse x_centers / y_centers computed above)
    xi, yi = np.meshgrid(x_centers, y_centers)

    # 4. Plotting
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Main 3D Surface
    # Transpose free_energy to match xi, yi dimensions
    surf = ax.plot_surface(xi, yi, free_energy.T, cmap='jet',
                           edgecolor='black', linewidth=0.1,
                           antialiased=False, rstride=1, cstride=1, alpha=1.0)

    # 5. The Projected Contour (The "floor" map seen in target images)
    # This places a 2D projection at the bottom of the Z-axis
    offset = 0 # Adjust if you want the floor higher or lower
    cset = ax.contourf(xi, yi, free_energy.T, zdir='z', offset=offset, cmap='jet', levels=20)

    # Labels and Scaling
    ax.set_xlabel('PC1 ($\AA$)', fontweight='bold', fontsize=12)
    ax.set_ylabel('PC2 ($\AA$)', fontweight='bold', fontsize=12)
    ax.set_zlabel('Gibbs Free Energy (kJ/mol)', fontweight='bold', fontsize=12)

    # Setting axis limits to match your data
    ax.set_zlim(0, energy_cap)

    # Colorbar
    cbar = fig.colorbar(surf, shrink=0.6, aspect=12, pad=0.1)
    cbar.set_label('Gibbs Free Energy (kJ/mol)', fontweight='bold')

    # View Angle to match target image style
    ax.view_init(elev=30, azim=45)

    # Clean up panes
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    output_name = f'Publication_FEL_{system_name}.png'
    plt.savefig(output_name, dpi=300, bbox_inches='tight')

    # ============================================================
    # SAVE FEL GRID
    # ============================================================

    print("Saving FEL grid...")

    summary_file = f"{system_name}_fel_summary.txt"
    grid_file = f"{system_name}_fel_grid.csv"

    with open(grid_file, "w") as f:

        f.write("PC1_Angstrom,PC2_Angstrom,DeltaG_kJmol\n")

        for i in range(len(x_centers)):
            for j in range(len(y_centers)):

                energy = analysis_energy[i, j]

                if np.isnan(energy):
                    continue

                f.write(
                    f"{x_centers[i]:.6f},"
                    f"{y_centers[j]:.6f},"
                    f"{energy:.6f}\n"
                )

    # ============================================================
    # SAVE SUMMARY
    # ============================================================

    print("Saving FEL summary...")

    with open(summary_file, "w") as f:

        f.write("FREE ENERGY LANDSCAPE SUMMARY\n")
        f.write("="*45+"\n\n")

        f.write(f"Frames analysed : {len(x)}\n")
        f.write(f"Histogram bins  : {nbins} x {nbins}\n")
        f.write(f"Temperature     : {T:.1f} K\n\n")

        f.write("GLOBAL MINIMUM\n")
        f.write("---------------------------\n")
        f.write(f"DeltaG : {global_min_energy:.4f} kJ/mol\n")
        f.write(f"PC1    : {global_pc1:.4f} Å\n")
        f.write(f"PC2    : {global_pc2:.4f} Å\n\n")

        f.write("FREE ENERGY STATISTICS\n")
        f.write("---------------------------\n")
        f.write(f"Mean ΔG : {mean_energy:.4f} kJ/mol\n")
        f.write(f"Max ΔG  : {max_energy:.4f} kJ/mol\n")
        f.write(f"SD ΔG   : {std_energy:.4f} kJ/mol\n\n")

        f.write("PC1 RANGE\n")
        f.write("---------------------------\n")
        f.write(f"Minimum : {pc1_min:.4f} Å\n")
        f.write(f"Maximum : {pc1_max:.4f} Å\n")
        f.write(f"Range   : {pc1_range:.4f} Å\n\n")

        f.write("PC2 RANGE\n")
        f.write("---------------------------\n")
        f.write(f"Minimum : {pc2_min:.4f} Å\n")
        f.write(f"Maximum : {pc2_max:.4f} Å\n")
        f.write(f"Range   : {pc2_range:.4f} Å\n")

    plt.close(fig)
    print(f"Final publication plot saved as: {output_name}")
    print(f"Saved: {summary_file}")
    print(f"Saved: {grid_file}")
    print("FEL analysis completed.")

if __name__ == "__main__":
    # Replace with your actual file names
    plot_publication_fel('pc1.xvg', 'pc2.xvg', 'C11')