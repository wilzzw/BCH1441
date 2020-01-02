# Graduate Unit Design

cd "C:/Users/Wilson/Documents/BCH1441"

# Test block #
set PDB_CODE 6MSM
set residue_name ALA

# Idea: read from JSON for each type of residue to decide what kind of atoms to show
set selection_atoms "N CA CB C"

# Load the structure file
#mol new ${path}trajectories/first_frame.gro type {gro} first 0 last -1 step 1 waitfor -1
mol new $PDB_CODE.pdb type {pdb} first 0 last -1 step 1 waitfor -1
mol modstyle 0 top NewCartoon
mol showrep top 0 off

# Select all such residues and use the alpha carbons to get the residue numbers
# This returns a list of residue numbers
# Also get the indices
set select_all_residues [atomselect top "protein and name CA and resname $residue_name"]
set all_resIDs [$select_all_residues get resid]
set all_residues [$select_all_residues get residue]

# The number of such residue
set num_residues [llength $all_residues]

# Create as many copies of the structure as there are such residue in the protein
# And show their residue representations as well as their surroundings
# And align them to a representative residue structure, which will be the first structure
for {set i 0} {$i < $num_residues} {incr i} {
    mol new $PDB_CODE.pdb type {pdb} first 0 last -1 step 1 waitfor -1
    set residue_to_show [lindex $all_residues $i]
    set resID_to_show [lindex $all_resIDs $i]
    mol modselect 0 top "protein and same residue as (within 4.5 of residue $residue_to_show)"
    mol addrep top
    mol modselect 1 top "protein and residue $residue_to_show"
    mol modstyle 1 top Licorice
    mol rename top ${residue_name}${resID_to_show}
    if {$i > 0} {
        set toalign_residue [atomselect top "protein and residue $residue_to_show and name $selection_atoms"]
        set transform [measure fit $toalign_residue $ref_residue]
        set full_protein [atomselect top "all"]
        $full_protein move $transform
        mol off top
    } else {
        set ref_residue [atomselect top "protein and residue $residue_to_show and name $selection_atoms"]
    }
}

# Use the first of such residue for the reference view
display resetview

# Measure dihedrals

# Cluster

# Set up lists for all clusters

# For each list, align them to a representative residue structure

# Combine and compute density. Refer to MDAnalysis source code or VMD source code for it.

# Show density