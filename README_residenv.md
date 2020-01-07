v1.0 completed 2020-01-07 by Zhi Wei (Wilson) Zeng
If you have any questions or problems with using the package, please contact: wilson.zeng@mail.utoronto.ca

residenv is a VMD script used for the intention of:
- visualizing residue environments in terms of densities of distributions of atom types
- intuitively visualizing residue chi-angle distributions in terms of densities of such distributions
- compare dihedral angles or surrounding environments of the same type of residues

Potential utilities:
- Rationalize protein residue functions in the context of its environments
- Observe dihedral angle biases
- Educational value of visualizing chi-angles
- Drug design

Installation & Description of package components
- To download and install VMD, go to: https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD (It should be free for academically affiliated)
- Example pdb files "6msm.pdb" and "5w81.pdb" are included
- "residenv.tcl": the main script of the program
- "tcllib" utilities: cloned from git@github.com:tcltk/tcllib.git. It contains parser for JSON files. See Lines 39-56 in "residenv.tcl" for more details
- "top_all27_prot_lipid_withHIS.inp": residue topology file required by the psfgen tcl-package to interpret residues and add hydrogens. This is the modified version of VMD "top_all27_prot_lipid.inp", which no longer recognizes HIS as a residue name, but instead one of HSD, HSE, and HSP, which will result in pdb reading error. Therefore, I copied the section that correspond to "HSD" and rename it as "HIS" and added it to the file, making "HSD" the default for "HIS". For more details, please inspect this INP file.
- "atomsAlign.json": stored PDB-names for the atoms of each residue to superimpose depending on viewing purposes.
- "atomToLookAt.json": stored PDB-names for the suggested atoms of each residue to visualize chi-angle distributions/densities
- "chiAngleAtoms.json": sotred PDB-names for atoms used to define atoms and selection in VMD for calculating chi-angles

How to use:
- Change the global variables in the script "residenv.tcl" will allow setting changes
- In VMD, in the command line interface OR in the "Tk Consoles" under "Extensions", simply type "source $path/residenv.tcl", where $path is the location to which this package is extracted to. This path will also be used as the variable HOME_PATH in "residenv.tcl". Please change it according to your situation.
- Within "residenv.tcl":
    > PDB_CODE is the pdb file you want to use as the dataset for analysis. Note: for VMD versions prior to 1.9.3, automatic downloads of pdb file based on the 4-character PDB code does not work. Please download the PDB file into $path.
    > RESIDUE_NAME is the name of the protein residue to be analyzed. Currently the accepted format is the three-letter amino acid code
    > MODE: If you want to visualize the chi-angle distributions, MODE will be "chi1", "chi2", ... etc. If you want to visualize the surrounding environment of the residue, MODE will be "env"
    > Basic filtering based on chi-angles is possible. For example, if we only want to consider residues with chi1 angles within a certain range, based on histograms, use the proc un-commented from Lines 170-179. Unfortunately, tcl does not have a good intuitive plotting method I can understand, so I recommend base your clustering decisions on chi-angle distribution plots created externally
    > Other default settings can be modified in script as well (e.g. resolution of density surface)
    > Don't make the script stop you from interacting/intervening with the VMD GUI!

Future development plans:
- Combine information from multiple PDB files
- More advanced filtering options. Currently there is only chi-angle filtering possible
    (Might be nice to be able to filter by secondary structures, B-factor, local resolution, occupancy, user-determined parameters, and even bioinformatic annotations)
- Include consideration of non-proteogenic chemical moieties & with reference to chemical fragment libraries
- Extend the package to Pymol and UCSF Chimera