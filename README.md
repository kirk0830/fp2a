# fp2a: First principle simulation softwares input script conversion to ABACUS
for input script conversion

### usage
input parameter description:
&basic_conversion
    mode = 'single', 'multiple' or 'examples'
            # single: convert a single qe input script
            # multiple: convert a set of qe input scripts
            # examples: convert a set of qe input scripts in the examples folder, in this case,
            #           the input scripts are organized in shell scripts and each one ends with
            #           'EOF'.
    script_dir = 'path to the folder containing the qe input script(s)'
    script_name = 'qe input script name'
            # if mode = 'single', then script_name is the name of the qe input script
            # if mode = 'multiple', then script_name is a cascade of qe input script names, separated by ','
            # if mode = 'examples', then this keyword is omitted and all files named as 'run_example' will be converted
    qe_version = 6.8
/# end of &basic_conversion

&overwrite_keywords
# it is very possible that you need to specify some special keywords for abacus
# for example, you may not want to use the same folder storing both pseudopotentials of qe and abacus
# so you can define pseudopotential of abacus here, if so, it will overwrite the one read from qe input script
# or the same folder will be used
    #pseudo_dir = './'
    #orbital_dir = './'
/# end of &overwrite_keywords

&additional_keywords
# well, if you use additional keywords, may be the comparison between qe and abacus
# is spoiled, so, make sure you know what you are doing.
    [abacus keywords] [values]
/# end of &additional_keywords
