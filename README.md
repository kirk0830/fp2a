# fp2a: First principle simulation softwares input script conversion to ABACUS
for input script conversion

## usage  
input parameter description:  
&basic_conversion  
    mode = 'single', 'multiple' or 'examples'  
            (comment) single: convert a single qe input script  
            (comment) multiple: convert a set of qe input scripts  
            (comment) examples: convert a set of qe input scripts in the examples folder, in this case,  
            (comment)           the input scripts are organized in shell scripts and each one ends with  
            (comment)           'EOF'.  
    script_dir = 'path to the folder containing the qe input script(s)'  
    script_name = 'qe input script name'  
            (comment) if mode = 'single', then script_name is the name of the qe input script  
            (comment) if mode = 'multiple', then script_name is a cascade of qe input script names, separated by ','  
            (comment) if mode = 'examples', then this keyword is omitted and all files named as 'run_example' will be converted  
    qe_version = 6.8  
/# end of &basic_conversion  
  
&overwrite_keywords  
(comment) it is very possible that you need to specify some special keywords for abacus  
(comment) for example, you may not want to use the same folder storing both pseudopotentials of qe and abacus  
(comment) so you can define pseudopotential of abacus here, if so, it will overwrite the one read from qe input script  
(comment) or the same folder will be used  
    (comment)pseudo_dir = './'  
    (comment)orbital_dir = './'  
/(comment) end of &overwrite_keywords  
  
&additional_keywords  
(comment) well, if you use additional keywords, may be the comparison between qe and abacus  
(comment) is spoiled, so, make sure you know what you are doing.  
    [abacus keywords] [values]  
/(comment) end of &additional_keywords  
