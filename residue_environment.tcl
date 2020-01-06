# Graduate Unit Design

# Global variables
set HOME_PATH "C:/Users/Wilson/Documents/BCH1441"
set PDB_CODE 6MSM
set RESIDUE_NAME LYS

set FILTER_LANGUAGE "protein and not resname UNK"
# To align GLY, MET, and LYS properly, hydrogens have to be added
set RLIST_NEEDH_TOALIGN {GLY MET LYS}

cd $HOME_PATH
# Make temporary folder
set TEMP_FOLDER .residue_env

# Damage control. If TEMP_FOLDER already exists...
if {[file isdirectory $TEMP_FOLDER]} {
    puts "Folder $TEMP_FOLDER already exists. Please change the name of TEMP_FOLDER"
} else {
    file mkdir $TEMP_FOLDER
}

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
set tclJsonPATH "${HOME_PATH}/tcllib/modules/json"
lappend auto_path $tclJsonPATH
package require json

# Function that reads JSON files
# Argument: jsonName: JSON file name with extension .json
proc readJSON {jsonName} {
    global HOME_PATH
    set fileContent [open "${HOME_PATH}/${jsonName}" r]
    set read_fileContent [read $fileContent]
    close $fileContent
    return [::json::json2dict $read_fileContent]
}

# Add hydrogens with the psfgen package
package require psfgen
proc addHydrogens {pdbFilename} {
    set topologyfile "top_all27_prot_lipid_withHIS.inp"
    # Interpret topology definitions
    topology $topologyfile
    set segid A
    segment $segid {
        first NONE
        last NONE
        pdb $pdbFilename
    }
    coordpdb $pdbFilename $segid
    guesscoord
    writepdb ${pdbFilename}
    resetpsf
}

# Load the structure file
# If addh is True, add hydrogens
proc strucLoad {addh} {
    global PDB_CODE
    global FILTER_LANGUAGE
    mol new ${PDB_CODE}.pdb type {pdb} first 0 last -1 step 1 waitfor -1
    if {addh} {
        puts "Note: hydrogens will be added to the pdb based on reasonable bond geometries"
        # Filter off residues that psfgen cannot read and save
        [atomselect top $FILTER_LANGUAGE] writepdb ${PDB_CODE}.pdb
        # Delete the non-hydrogenated one
        mol delete top
        # Add hydrogen and overwrite
        puts "Overwriting ${PDB_CODE}.pdb..."
        addHydrogens ${PDB_CODE}.pdb
        # Re-load the structure file
        mol new ${PDB_CODE}.pdb type {pdb} first 0 last -1 step 1 waitfor -1
    }
}

# Read "atomsAlign.json"
set dict_atomsAlign [readJSON "atomsAlign.json"]

# Get atoms selected for superposition for the residue type
#set selection_atoms "N CA CB C"
set selection_atoms [dict get $dict_atomsAlign $RESIDUE_NAME]

# Load the structure file
# lsearch returns the index of the occurance of an element in a list (0-based)
# If the element does not exist in the list, it returns -1
# But in tcl, greater than zero is True... So +1.
set needAddHydrogen [lsearch $RLIST_NEEDH_TOALIGN $RESIDUE_NAME]
# If the residue is GLY, LYS, or MET, add hydrogens to the structure
strucLoad [expr $needAddHydrogen + 1]

mol modstyle 0 top NewCartoon
mol showrep top 0 off

# Select all such residues and use the alpha carbons to get the residue numbers
# This returns a list of residue numbers
# Also get the indices
### Will have other treatments if there are dihedral criteria
set select_all_residues [atomselect top "protein and name CA and resname $RESIDUE_NAME"]
set all_resIDs [$select_all_residues get resid]
set all_residues [$select_all_residues get residue]

# How many of such residues are there?
set num_residues [llength $all_residues]

# Create as many copies of the structure as there are such residue in the protein
# And show their residue representations as well as their surroundings
# And align them to a representative residue structure, which will be the first structure
for {set i 0} {$i < $num_residues} {incr i} {
    # Load copies of the PDB structure and rename them to all residue names
    mol new ${PDB_CODE}.pdb type {pdb} first 0 last -1 step 1 waitfor -1

    set resindex_to_show [lindex $all_residues $i]
    set resID_to_show [lindex $all_resIDs $i]
    mol rename top ${RESIDUE_NAME}${resID_to_show}

    # Showing subsets
    # Show the residue and its surrounding residues (within 4.5 Ångstrom of heavy atom cut-off) in wire/line representation (default)
    # Also show the selected residue itself in licorice representation
    # Selection language for residue and its surroundings
    # Future development includes non-protein ligands. Require necessary respective .inp files
#    set selang_contact "protein and same residue as (within 4.5 of residue $resindex_to_show)"
# Added to evade error with UNK, which is proteogenic with unknown identity (modelled as poly-Ala)
# This is also why currently I cannot read the whole .pdb file with psfgen and add hydrogens for once and for all..
# Consider writing as a filter function
    set selang_contact "protein and not resname UNK and same residue as (within 4.5 of residue $resindex_to_show)"
    mol modselect 0 top $selang_contact
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

    # Save individual residue and its surrounding into PDB for later merging and computing densities
    # Separately, so that we can distinguish GLU and surrounding GLU during density calculation
    [atomselect top "$selang_contact and not residue $resindex_to_show"] writepdb ${TEMP_FOLDER}/${RESIDUE_NAME}${resID_to_show}_ENV.pdb
    [atomselect top "residue $resindex_to_show"] writepdb ${TEMP_FOLDER}/${RESIDUE_NAME}${resID_to_show}.pdb

}

# Reference view. Get close to the view from the perspective of the residue
display resetview

# Measure dihedrals

# Chi angle clustering recommendations

# Cluster via histogram

# Set up lists for all clusters

proc processPDB_writeout {pdbfile} {
    set pdbOpen [open $pdbfile r]
    set pdbContent [read $pdbOpen]
    close $pdbOpen
    set pdbLines [split $pdbContent "\n"]

    set numberOfLines [llength $pdbLines]
    set start 1
    set end [expr $numberOfLines - 3]

    return [lrange $pdbLines $start $end]
}


set mergedEnv [open "merged_${RESIDUE_NAME}_ENV.pdb" w]
puts $mergedEnv "TITLE     FUSED PDB FOR DENSITY CALCULATION"
foreach pdbfile [glob ${TEMP_FOLDER}/*_ENV.pdb] {
    set lines2write [processPDB_writeout $pdbfile]
    foreach line $lines2write {
        puts $mergedEnv $line
    }
    file delete $pdbfile
}
puts $mergedEnv "END"
close $mergedEnv

set mergedResidue [open "merged_${RESIDUE_NAME}.pdb" w]
puts $mergedResidue "TITLE     FUSED PDB FOR DENSITY CALCULATION"
foreach pdbfile [glob ${TEMP_FOLDER}/*.pdb] {
    set lines2write [processPDB_writeout $pdbfile]
    foreach line $lines2write {
        puts $mergedResidue $line
    }
    file delete $pdbfile
}
puts $mergedResidue "END"
close $mergedResidue

file delete $TEMP_FOLDER

# Compute and show density
# Something like...
#volmap density [atomselect top "oxygen"] -res 0.5 -weight mass