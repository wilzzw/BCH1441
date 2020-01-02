# Graduate Unit Design

# Test block #
set HOME_PATH "C:/Users/Wilson/Documents/BCH1441"
set PDB_CODE 6MSM
set RESIDUE_NAME GLY
set tclJsonPATH "${HOME_PATH}/tcllib/modules/json"

cd $HOME_PATH

# Idea: read from JSON for each type of residue to decide what kind of atoms to show
# Select a subset of atoms (typically four non-coplanar atoms) for aligning the residues
# This will depend on the residue type specified in $RESIDUE_NAME. Specific atom names to be used for alignment is stored in JSON file "atomsAlign.json"
# Here a package is used to read JSON. It appears that unless I installed tcl/tk properly (not as an extension of VMD), it won't find proper package to read JSON
# So the following steps are needed:

### Step 1: Make a local repository clone of git@github.com:tcltk/tcllib.git (getting the source codes of the package)
### -- e.g. in Linux/Unix (mine is Bash on Ubuntu on Windows):
### -- mkdir tcllib; cd tcllib;
### -- git init; git remote add origin git@github.com:tcltk/tcllib.git;
### -- git pull origin master
###
### Step 2: Add the absolute path to "tcllib/modules/json" to the tcl environment variable $auto_path in VMD
### -- For me, the command line I type in the tcl/tk console of VMD is:
### -- lappend auto_path "C:/Users/Wilson/Documents/BCH1441/tcllib/modules/json"
###
### Step 3: Loading the package json in tcl should work now:
### -- package require json

# Step 2 and 3 have already been implemented in the code right underneath. User has to do Step 1 by themselves and change the path to "json" package accordingly
# I put that in the global variable $tclJsonPATH
lappend auto_path $tclJsonPATH
package require json

# Read "atomsAlign.json"
set atomsAlign [open "${HOME_PATH}/atomsAlign.json" r]
set read_atomsAlign [read $atomsAlign]
close $atomsAlign
set dict_atomsAlign [::json::json2dict $read_atomsAlign]

# Get atoms selected for superposition for the residue type
#set selection_atoms "N CA CB C"
set selection_atoms [dict get $dict_atomsAlign $RESIDUE_NAME]

# Load the structure file
#mol new ${path}trajectories/first_frame.gro type {gro} first 0 last -1 step 1 waitfor -1
mol new $PDB_CODE.pdb type {pdb} first 0 last -1 step 1 waitfor -1
mol modstyle 0 top NewCartoon
mol showrep top 0 off

# Select all such residues and use the alpha carbons to get the residue numbers
# This returns a list of residue numbers
# Also get the indices
### Will have other treatments if there are dihedral criteria
set select_all_residues [atomselect top "protein and name CA and resname $RESIDUE_NAME"]
set all_resIDs [$select_all_residues get resid]
set all_residues [$select_all_residues get residue]

# The number of such residue
set num_residues [llength $all_residues]

# Create as many copies of the structure as there are such residue in the protein
# And show their residue representations as well as their surroundings
# And align them to a representative residue structure, which will be the first structure
for {set i 0} {$i < $num_residues} {incr i} {
    # Load copies of the PDB structure and rename them to all residue names
    mol new $PDB_CODE.pdb type {pdb} first 0 last -1 step 1 waitfor -1
    set resindex_to_show [lindex $all_residues $i]
    set resID_to_show [lindex $all_resIDs $i]
    mol rename top ${RESIDUE_NAME}${resID_to_show}

    # Showing subsets
    # Show the residue and its surrounding residues (within 4.5 Ångstrom of heavy atom cut-off) in wire/line representation (default)
    # Also show the selected residue itself in licorice representation
    mol modselect 0 top "protein and same residue as (within 4.5 of residue $resindex_to_show)"
    mol addrep top
    mol modselect 1 top "protein and residue $resindex_to_show"
    mol modstyle 1 top Licorice
    
    # The first of such residues will be the reference to align the rest of the residues to
    # Except for the first of such residues, do not display the visualization (can be manually turned back on)
    if {$i > 0} {
        set toalign_residue [atomselect top "protein and residue $resindex_to_show and name $selection_atoms"]
        set transform [measure fit $toalign_residue $ref_residue]
        set full_protein [atomselect top "all"]
        $full_protein move $transform
        mol off top
    } else {
        set ref_residue [atomselect top "protein and residue $resindex_to_show and name $selection_atoms"]
    }
}

# Reference view. Get close to the view from the perspective of the residue
display resetview

# Measure dihedrals

# Cluster

# Set up lists for all clusters

# For each list, align them to a representative residue structure

# Combine and compute density. Refer to MDAnalysis source code or VMD source code for it.

# Show density