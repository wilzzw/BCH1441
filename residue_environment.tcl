# Graduate Unit Design

# Global variables
set HOME_PATH "C:/Users/Wilson/Documents/BCH1441"
set PDB_CODE 6MSM
set RESIDUE_NAME LYS
set MODE env

# Added to evade error with UNK, which is proteogenic with unknown identity (modelled as poly-Ala)
# This is also why I have to filter them out so that psfgen package can properly read the protein.
# Also, treating them as Alanines is pointless, misleading, and non-biological
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

    if {[file exists ${PDB_CODE}.pdb]} {
        mol new ${PDB_CODE}.pdb type {pdb} first 0 last -1 step 1 waitfor -1
    } else {
        mol pdbload $PDB_CODE
    }
    
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

proc selectAtoms {} {
    global RESIDUE_NAME
    global MODE
    # Read "atomsAlign.json"
    set dict_atomsAlign [readJSON "atomsAlign.json"]

    # Get atoms selected for superposition for the residue type
    # This retrieves a dictionary, whose keys are the "modes" of viewing
    # (e.g. viewing chi angle distributions/densities or looking at the residue's surrounding environment densities)
    set selection_atoms_allmodes [dict get $dict_atomsAlign $RESIDUE_NAME]
    return [dict get $selection_atoms_allmodes $MODE]
}

proc needAddHydrogens {atomNameList} {
    foreach aname $atomNameList {
        if {[string index $aname 0] == "H"} {
            return 1
        }
    }
    return 0
}

# Select atoms for residue structural superposition
set selected_atoms [selectAtoms]

# Load the structure file
# Add hydrogens if needed (i.e. when selected_atoms include hydrogens)
strucLoad [needAddHydrogens $selected_atoms]

# Set the whole protein in Cartoon representation and hide it
mol modstyle 0 top NewCartoon
mol showrep top 0 off

# Select all such residues and use the alpha carbons to get the residue numbers
# This returns a list of residue numbers
# Also get the indices
### Will have other treatments if there are dihedral criteria
set select_all_residues [atomselect top "protein and name CA and resname $RESIDUE_NAME"]
#set all_resIDs [$select_all_residues get resid]
set all_residues [$select_all_residues get residue]

# Measure dihedrals
proc measureDihedral {alist mol} {
    return [measure dihed $alist molid $mol]
}

# Measure chi-angles
proc measureChi {n resindex mol} {
    set atomsToDefine [lrange $atomsDefs [expr $n - 1] [expr $n + 3]]
    set atomList [list]
    foreach a $atomsToDefine {
        set possibleAtomNames [join $a " "]
        set atom [atomselect $mol "residue $resindex and name $possibleAtomNames"]
        set atomindex [$atom get index]
        lappend atomList $atomindex
    }
    return [measureDihedral $atomList $mol]
}

# Dummy condition for testing
#proc condition {resindex} {
 #   return 1
#}

# Chi1 filtering
proc condition {resindex} {
    set diheadralAngle [measureChi 1 $resindex top]
    if {($dihedralAngle >= 50) && ($dihedralAngle <= 70)} {
        return 1
    } else {
        return 0
    }
}

proc filterResidues {residueList} {
    set filteredList [list]
    foreach i $residueList {
        if {[condition $i]} {
            lappend filteredList $i 
        }
    }
    return $filteredList
}

set qualifiedResidues [filterResidues $all_residues]
set filteredSelection [atomselect top "name CA and residue $qualifiedResidues"]
set qualifiedResIDs [$filteredSelection get resid]

# How many of such residues are there?
set num_residues [llength $qualifiedResidues]

# Create as many copies of the structure as there are such residue in the protein
# And show their residue representations as well as their surroundings
# And align them to a representative residue structure, which will be the first structure
set heavy_cutoff 4.5
for {set i 0} {$i < $num_residues} {incr i} {
    # Load copies of the PDB structure and rename them to all residue names
    mol new $PDB_CODE.pdb type {pdb} first 0 last -1 step 1 waitfor -1

    set resindex_to_show [lindex $qualifiedResidues $i]
    set resID_to_show [lindex $qualifiedResIDs $i]
    mol rename top ${RESIDUE_NAME}${resID_to_show}

    # Showing subsets
    # Show the residue and its surrounding residues (heavy atom within 4.5 Ångstrom of heavy atom cut-off) in wire/line representation (default)
    # Also show the selected residue itself in licorice representation (more emphasized)
    # Selection language for residue and its surroundings. No non-protein moieties are currently involved
    # Future development includes non-protein ligands. Require necessary respective .inp files
    set selang_contact "same residue as (mass > 1 and within $heavy_cutoff of residue $resindex_to_show)"
    mol modselect 0 top $selang_contact
    mol addrep top
    mol modselect 1 top "protein and residue $resindex_to_show"
    mol modstyle 1 top Licorice

    # The first of such residues will be the reference to superimpose the rest of the residues to
    # Except for the first of such residues, do not display the visualization (can be manually turned back on through GUI)
    if {$i == 0} {
        set ref_residue [atomselect top "residue $resindex_to_show and name $selected_atoms"]
    } else {
        set toalign_residue [atomselect top "residue $resindex_to_show and name $selected_atoms"]
        set transform [measure fit $toalign_residue $ref_residue]
        set full_protein [atomselect top "all"]
        $full_protein move $transform
        mol off top
    }

    # Save individual residue and its surrounding into PDB for later merging and computing densities
    # Save the residue itself and its surroundings separately, so that we can distinguish the residue and its surrounding during density calculation
    [atomselect top "$selang_contact and not residue $resindex_to_show"] writepdb ${TEMP_FOLDER}/${RESIDUE_NAME}${resID_to_show}_ENV.pdb
    [atomselect top "residue $resindex_to_show"] writepdb ${TEMP_FOLDER}/${RESIDUE_NAME}${resID_to_show}.pdb

}

# Reference view. Get close to the view from the perspective of the residue
display resetview

# Chi angle clustering recommendations
# Cluster via histogram using external programs

# Set up lists for all clusters/ filtering function

# Read lines from a PDB file saved by VMD
# Returns the second line till the line before the "END" line
# Used as part of the utility to write several PDB files into a merged PDB file
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

proc mergePDB {suffix} {
    global RESIDUE_NAME
    global TEMP_FOLDER

    if {[string length $suffix] > 0} {
        set suffix "_$suffix"
    }

    set merged [open "merged_${RESIDUE_NAME}${suffix}.pdb" w]
    puts $merged "TITLE     FUSED PDB FOR DENSITY CALCULATION"

    foreach pdbfile [glob ${TEMP_FOLDER}/*${suffix}.pdb] {
        set lines2write [processPDB_writeout $pdbfile]
        foreach line $lines2write {
            puts $merged $line
        }
        file delete $pdbfile
    }
    puts $merged "END"
    close $merged
    return "merged_${RESIDUE_NAME}${suffix}.pdb"
}

# Write merged PDB for the surroundings of requested residues
set mergedEnv [mergePDB ENV]
# Write merged PDB for the requested residues themselves
set mergedResidues [mergePDB]

file delete $TEMP_FOLDER

# Load the merged environment if the MODE is env
# Otherwise, it makes more sense to load and analyze the residue itself
if {$MODE == env} {
    mol new $mergedEnv type {pdb} first 0 last -1 step 1 waitfor -1
} else {
    mol new $mergedResidues type {pdb} first 0 last -1 step 1 waitfor -1
}

#file delete $mergedEnv
#file delete $mergedResidues

# Compute and show density
proc calcDensity {selection resolution outputName} {
    volmap density [atomselect top $selection] -res $resolution -weight mass -o $outputName
    mol new $outputName type {dx} first 0 last -1 step 1 waitfor -1
}

proc atomToLookAt {} {
    global MODE
    if {$MODE == "env"} {
        return
    }
    set atomsOfInterest [readJSON "atomToLookAt.json"]
    set atomsSelected [dict get [dict get $atomsOfInterest $RESIDUE_NAME] $MODE]
    puts "Recommended atoms to compute densities for: $atomsSelected"
    set atomsSelectedLang [join $atomsSelected " "]
    return "name $atomsSelectedLang"
}

set candidateAtoms [atomToLookAt]
#calcDensity "oxygen" 0.5 test.dx
calcDensity $candidateAtoms 0.5 test.dx