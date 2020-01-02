# Graduate Unit Design

cd "C:/Users/Wilson/Documents/BCH1441"

# Test block #
set PDB_CODE 6MSM
set residue PHE

# Load the structure file
#mol new ${path}trajectories/first_frame.gro type {gro} first 0 last -1 step 1 waitfor -1
mol new $PDB_CODE.pdb type {pdb} first 0 last -1 step 1 waitfor -1

# Select all such residues
set all_residues [atomselect 0 "protein and resname $residue"]

# Measure dihedrals


# Cluster

# Set up lists for all clusters

# For each list, align them to a representative residue structure
# Load multiple copies, or just copy from top
# For each of them, pick out only the residues within 4.5 Angstroms and show them

# Combine and compute density. Refer to MDAnalysis source code or VMD source code for it.

# Show density