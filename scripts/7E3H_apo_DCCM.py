import os
import numpy as np
import MDAnalysis as mda


def calculate_dccm_concatenated():

    # ============================================================
    # Input files
    # ============================================================
    topology_file = "topology_fix.pdb"
    trajectory_file = "concat_fit.xtc"

    # Output files
    matrix_csv = "dccm_matrix.csv"
    residue_txt = "residue_ids.txt"

    # ============================================================
    # Check input files
    # ============================================================
    if not os.path.exists(topology_file):
        print(f"Error: {topology_file} not found.")
        return

    if not os.path.exists(trajectory_file):
        print(f"Error: {trajectory_file} not found.")
        return

    # ============================================================
    # Load trajectory
    # ============================================================
    print(f"Loading Concatenated System: {trajectory_file}...")
    u = mda.Universe(topology_file, trajectory_file)

    # ============================================================
    # Select C-alpha atoms
    # ============================================================
    ca = u.select_atoms("name CA")

    print(f"Found {len(ca)} C-alpha atoms.")

    # ============================================================
    # Extract coordinates
    # ============================================================
    print("Extracting coordinates...")

    coordinates = np.array(
        [ca.positions.copy() for ts in u.trajectory[::10]],
        dtype=np.float32
    )

    n_frames = coordinates.shape[0]
    n_atoms = coordinates.shape[1]

    print(f"Frames used: {n_frames}")
    print(f"Residues: {n_atoms}")

    # ============================================================
    # Calculate displacements
    # ============================================================
    print("Calculating displacements...")

    mean_pos = coordinates.mean(axis=0)
    disp = coordinates - mean_pos

    disp = disp.reshape(n_frames, n_atoms * 3)

    # ============================================================
    # Covariance matrix
    # ============================================================
    print("Calculating covariance matrix...")

    cov = np.dot(disp.T, disp) / n_frames

    # ============================================================
    # Calculate DCCM
    # ============================================================
    print("Calculating DCCM...")

    cor_mat = np.zeros((n_atoms, n_atoms), dtype=np.float32)

    for i in range(n_atoms):

        if i % 25 == 0:
            print(f"Residue {i+1}/{n_atoms}")

        dot_ii = np.trace(cov[i*3:i*3+3, i*3:i*3+3])

        for j in range(n_atoms):

            dot_jj = np.trace(cov[j*3:j*3+3, j*3:j*3+3])
            dot_ij = np.trace(cov[i*3:i*3+3, j*3:j*3+3])

            denominator = np.sqrt(dot_ii * dot_jj)

            if denominator > 0:
                cor_mat[i, j] = dot_ij / denominator
            else:
                cor_mat[i, j] = 0.0

    # ============================================================
    # Save raw data
    # ============================================================
    print("Saving raw data...")

    np.savetxt(
        matrix_csv,
        cor_mat,
        delimiter=",",
        fmt="%.6f"
    )

    np.savetxt(
        residue_txt,
        ca.resids,
        fmt="%d"
    )

    print(f"Saved: {matrix_csv}")
    print(f"Saved: {residue_txt}")

    print("\n======================================")
    print("DCCM calculation completed successfully")
    print("======================================")
    print("Generated files:")
    print("  dccm_matrix.csv")
    print("  residue_ids.txt")


if __name__ == "__main__":
    calculate_dccm_concatenated()