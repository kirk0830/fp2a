# FP2A: Convert input script of first-principle software to ABACUS
Author: Yike Huang, member of ABACUS developer team.
## Basic usage  
Just need to enter one-line command in your prompt:  
`python fp2a_driver.py`  
You would be better store all of them in the same directory, although storing them seperately will certainly work.
## Utilization context
This package is basically for generating ABACUS input scripts from converting examples of Quantum ESPRESSO. The conversion in batch has been tested for PW/examples/example01..14. But there are still keywords failed to convert to the one corresponding in ABACUS. Presently the RISM (solvation module of QE), OCCUPATION (to force occupation of atoms) are not implemented for conversion.
## Detailed parameter description  
### About input script
The input file of fp2a package is set as fp2a.inp and is not expected to change (although it can be changed in fp2a_io.py). The contents of fp2a.inp is organized in following way:  
#### Section
`basic_conversion`
#### Description
This section is for defining basic information of conversion.
#### Contents
    &basic_conversion  
    mode = 'single'
Firstly it is needed to define the mode of conversion, if only 1 file is needed to convert, or multiple files, or need to convert examples provided by Quantum ESPRESSO. Therefore avaiable choices are: `single`, `multiple` and `examples`.  

    script_dir =   
    script_name =  
Then define the path and name of QE input script. For not `single` mode, e.g. the `multiple` mode, you can define not only one file name for keyword `script_name`. To seperate multiple names of file, use comma `,`. For `example` mode `script_name` keyword is omitted automatically.

    qe_version = 6.8  
Due to there is a big change in the input script organization since the publish of QE 7.1, DFT+U module now is called by setting an individual `HUBBARD` section instead of setting several `Hubbard_U(i)`. Therefore it is needed to identify different versions of QE.  

    /   
a `/` ends a section. Then it is the end of present `basic_conversion` section, and turn to another.   
#### Section
`overwrite_keywords`
#### Description
Define ABACUS keywords and its value to overwrite.  
#### Contents

    &overwrite_keywords   
    pseudo_dir = './'  
    orbital_dir = './'  
    / 
Similarly it is very possible for users expecting to adjust the real directory for running ABACUS, storing generated files. In ABACUS, basically this is controlled by two keywords above. If define any keywords in this section, corresponding keywords will be overwrite.   
#### Section
`additional_keywords`
#### Description
Define ABACUS keywords and append to converted input script.  
#### Contents

    &additional_keywords  
    ...
In this  section it is possible to switch from planewave basis function (`basis_type pw`) to lcao one (`basis_type lcao`), i.e., the atomic numerical orbital, for achieving low scaling of $O(N)$. To switch to `lcao`, use `numerical_orbitals` section:

        &numerical_orbitals
            use_nao = .false.
The keyword `use_nao` determines whether use "nao" (Numerical atomic orbital) to expand Kohn-Sham orbitals. If set it to `.true.`, program will switch the system from `basis_type pw` to `basis_type lcao` automatically and, add corresponding numerical orbitals information into `STRU` the ABACUS structure file.  
___Important!___ Next one needs to define the path where numerical orbitals are stored. For downloading the lastest full set (as long as it is available) of numerical orbitals, visit our official website http://abacus.ustc.edu.cn/pseudo/list.htm. But one version of AllOrbitals 2.0 is also provided in this repository. Its nomenclature composes information of both pseudopotential and orbitals. Say "SG15-Version1p0" is the "SG15"-type pseudopotential with its version 1.0. For "AllOrbitals-Version2p0" it is the second edition of AllOrbitals.

            nao_dir = ./SG15-Version1p0__AllOrbitals-Version2p0/
For atomic orbitals it is also important to declare its composition. ABACUS provides two kinds of valence-split basis function, one is `dzp` and the other is `tzdp`.

            basis_type = dzp
Even for one single element and a given functional (although in AllOrbitals-Version2p0, only `gga` is supported), there are not only one kind of numerical orbitals, it is because different cutoff can be chosen and will affect accurancy and efficiency of calculation. In ABACUS, numerical orbitals differ mainly in cutoff of radius instead of energy. So use `radius_min` or `radius_max` is always a good choice for keyword `select_nao`, although selection of numerical orbitals by energy cutoff and even combination of these two are supported as well.

            select_nao = radius_min
(All available choices for keyword `select_nao` are: `radius_min`, `radius_max`, `energy_min`, `energy_max`, `energy_min_radius_min`, `energy_min_radius_max`, `energy_max_radius_min` and `energy_max_radius_max`.)

            atom_species = Pt, O, Ti

Define atomic species. If not specified here, program will read from QE input script the `ATOMIC_SPECIES` section.


            cutoff_list_radius = 0, 0, 0
            cutoff_list_energy = 100, 100, 100
After define keyword `select_nao`, it is, compulsory to know what do those available options mean.  
`radius_min`: program will select a numerical orbital with the minimal cutoff radius. However, the range for searching numerical orbital with the minimal radius cutoff should be confined, it is done by setting a minimal value of radius value. In practice if the threshold is set to 6 and originally there are available choices like 6, 7, 8, 9, 10, program will use \*\_gga\_ __6__ \_au\_100Ry\_\*.orb the numerical orbital file. If set threshold to 8, with `select_nao radius_min`, program will use \*\_gga_ __8__ _au\_100Ry\_\*.orb.

        /
    /
Do not forget end both `numerical_orbitals` and `additional_keywords` section.