'''
fp2a: First-principles to ABACUS
description:
this driver controls the flow of conversion from qe input script to the abacus one
an EXTERNAL INPUT SCRIPT named fp2a.inp is needed to control the flow of the conversion

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

'''
import fp2a_io
import fp2a_jsonRotation
import fp2a_qe2qeJson
import fp2a_qeJson2abacus
import os
import subprocess
import json
import time
def main():
    
    if not os.path.isfile('fp2a.inp'):
        print('fp2a.inp is not found, please check')
        exit()
    else:
        basicConversion, overwriteKeywords, additionalKeywords = fp2a_io.readInputScript()
        database_dir = basicConversion['database_dir']
        boolQeVersion_7_1 = False
        if float(basicConversion['qe_version']) >= 7.1:
            boolQeVersion_7_1 = True
        else:
            boolQeVersion_7_1 = False
        if not os.path.isfile(database_dir+'depuis_le_abacus.json'):
            print('depuis_le_abacus.json is not found, please check')
            exit()
        else:
            jsonDatabaseFileName = fp2a_jsonRotation.jsonRotate('qac', database_dir+'depuis_le_abacus.json')
            with open(jsonDatabaseFileName) as jsonDatabaseFile:
                jsonDatabase = json.load(jsonDatabaseFile)
            if basicConversion['mode'] == 'single':
                [qeJson] = fp2a_qe2qeJson.qe2qeJson(
                    basicConversion['script_dir']+basicConversion['script_name'][0], 
                    boolVersion7_1 = boolQeVersion_7_1,
                    parseMode = 'single',
                    boolDebugMode = False
                    )
                fp2a_qeJson2abacus.generateABACUSFiles(
                    mode = 'on-the-fly', qeInputScriptJson = qeJson,
                    keywordConversionJson = jsonDatabase, qeVersion_7_1 = boolQeVersion_7_1,
                    overwriteKeywords = overwriteKeywords, additionalKeywords = additionalKeywords
                )
            elif basicConversion['mode'] == 'multiple':
                for scriptName in basicConversion['script_name']:
                    [qeJson] = fp2a_qe2qeJson.qe2qeJson(
                        basicConversion['script_dir']+scriptName, 
                        boolVersion7_1 = boolQeVersion_7_1,
                        parseMode = 'single',
                        boolDebugMode = False
                        )
                    fp2a_qeJson2abacus.generateABACUSFiles(
                        mode = 'on-the-fly', qeInputScriptJson = qeJson, 
                        keywordConversionJson = jsonDatabase, qeVersion_7_1 = boolQeVersion_7_1,
                        overwriteKeywords = overwriteKeywords, additionalKeywords = additionalKeywords
                    )
                    subprocess.run('mv '+'INPUT'+' '+'INPUT-'+scriptName)
                    subprocess.run('mv '+'STRU'+' '+'STRU-'+scriptName)
                    subprocess.run('mv '+'KPT'+' '+'KPT-'+scriptName)
                    if os.path.isfile('KLINE'):
                        subprocess.run('mv '+'KLINE'+' '+'KLINE'+scriptName)
            elif basicConversion['mode'] == 'examples':
                print('For mode = examples, definition of input script is omitted')
                presentDirectory = os.getcwd()
                print("Validity check for input script parsing")
                print("Current directory: "+presentDirectory)
                directories = next(os.walk(presentDirectory))[1]

                for directory in directories:
                    print("Test runs in directory: "+directory)
                    files = [f for f in os.listdir(directory) if not f.startswith('.')]
                    for file in files:
                        if file.startswith("run_example"):
                            print("Test run for file: "+file)
                            if directory == "cluster_example": continue

                            qeJsons = fp2a_qe2qeJson.qe2qeJson(
                                f"{directory}\\run_example",
                                boolVersion7_1 = boolQeVersion_7_1,
                                parseMode = 'examples',
                                boolDebugMode = False
                                )
                            
                            # to parse file like input script defined in the shell script, must set mode = 'examples'
                            for qeJson in qeJsons:
                                fp2a_qeJson2abacus.generateABACUSFiles(
                                    mode = 'on-the-fly', qeInputScriptJson = qeJson,
                                    keywordConversionJson = jsonDatabase, qeVersion_7_1 = boolQeVersion_7_1,
                                    overwriteKeywords = overwriteKeywords, additionalKeywords = additionalKeywords
                                )
                                time.sleep(1)
                                if 'prefix' in qeJson['control']:
                                    identifier = qeJson['control']['prefix']
                                else:
                                    # use time stamp as identifier
                                    identifier = time.strftime("%Y%m%d-%H%M%S")
                                subprocess.run('mv '+'INPUT'+' '+'INPUT-'+identifier)
                                subprocess.run('mv '+'STRU'+' '+'STRU-'+identifier)
                                subprocess.run('mv '+'KPT'+' '+'KPT-'+identifier)
                                if os.path.isfile('KLINE'):
                                    subprocess.run('mv '+'KLINE'+' '+'KLINE-'+identifier)   
            else:
                print('mode is not defined correctly, please check')
                exit()
if __name__ == '__main__':
    main()
    print('--fp2a finished.')