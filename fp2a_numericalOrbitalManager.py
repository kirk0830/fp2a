import os
import re
import json
from fp2a_keywordsProcessing import keywordsWrite, keywordsRead

NAO_DATABASE = 'SG15-Version1p0__AllOrbitals-Version2p0'
NAO_JSON = 'fp2a_NAO.json'

def NAODatabaseGenerator(
    allOrbitalsDatabasePath = NAO_DATABASE, 
    NAONamePattern = r'[A-Z][a-z]?\w?_[a-z]{3}_([0-9]+)([\.]*)([0-9]*)au_([0-9]*)Ry_([0-9][a-z])+\.orb',
    boolSaveToJsonFile = False,
    jsonFileNameWithPath = NAO_JSON,
    debugMode = False
    ):
    '''
    This function is used to generate a dictionary of ABACUS NAOs whose name in format of
    [element]_[functional]_[radiusCutoff]au_[energyCutoff]Ry_[basisFunctions].orb\n
    The dictionary is saved in a json file named fp2a_NAO.json in the current directory by default.\n
    The dictionary is in the format of:\n
    NAOdictionary = {
        element: {
            functional: {
                
                    }
                }
            }
        }\n
    @param allOrbitalsDatabasePath: the path of the directory where the NAOs are stored\n
    @param NAONamePattern: the regular expression pattern of the NAO name\n
    @param boolSaveToJsonFile: whether to save the dictionary in a json file\n
    @param jsonFileNameWithPath: the path of the json file\n
    @param debugMode: whether to print the debug information\n
    @return: the dictionary of NAOs
    '''
    NAOdictionary = {}
    os.chdir(allOrbitalsDatabasePath)
    folderNames = [folder for folder in os.listdir() if os.path.isdir(folder)]
    for folderName in folderNames:
        words = re.split(pattern = "_", string = folderName)
        #atomic_index = words[0]
        element = words[1]
        zetaSelection = words[2].lower()
        os.chdir(folderName)
        fileNames = [file for file in os.listdir() if os.path.isfile(file)]
        for fileName in fileNames:
            match = re.search(NAONamePattern, fileName)
            if match:
                words = re.split(pattern = "_|.orb|au|Ry", string = fileName)
                element = keywordsRead(words[0])
                functional = keywordsRead(words[1])
                radiusCutoff = keywordsRead(words[2].replace('au', ''))
                energyCutoff = keywordsRead(words[4].replace('Ry', ''))
                basisFunctions = keywordsRead(words[6])
                if debugMode:
                    print("fileName: ", keywordsWrite(fileName))
                    print("element: ", keywordsWrite(element))
                    print("functional: ", keywordsWrite(functional))
                    print("radiusCutoff: ", keywordsWrite(radiusCutoff))
                    print("energyCutoff: ", keywordsWrite(energyCutoff))
                    print("basisFunctions: ", keywordsWrite(basisFunctions))
                if element not in NAOdictionary:
                    NAOdictionary[element] = {}
                if functional not in NAOdictionary[element]:
                    NAOdictionary[element][functional] = {}
                if zetaSelection not in NAOdictionary[element][functional]:
                    NAOdictionary[element][functional][zetaSelection] = {}
                if fileName not in NAOdictionary[element][functional][zetaSelection]:
                    NAOdictionary[element][functional][zetaSelection][fileName] = {}
                NAOdictionary[element][functional][zetaSelection][fileName] = {
                    "radiusCutoff": radiusCutoff,
                    "energyCutoff": energyCutoff,
                    "basisFunctions": basisFunctions
                    }
        os.chdir('..')
    os.chdir('..')
    if boolSaveToJsonFile:
        with open(jsonFileNameWithPath, 'w') as fp:
            json.dump(NAOdictionary, fp, indent = 4)
    return NAOdictionary

def openNAODatabase():
    if os.path.isfile(NAO_JSON):
        with open(NAO_JSON, 'r') as fp:
            NAOdictionary = json.load(fp)
        return NAOdictionary
    else:
        return NAODatabaseGenerator(
            boolSaveToJsonFile = True
        )

'''
def findNAOByElement(element):
    NAOdictionary = openNAODatabase()
    if element in NAOdictionary:
        returnList = []
        for functional in NAOdictionary[element]:
            returnList += list(NAOdictionary[element][functional][zetaSelection].keys())
        return returnList
    else:
        print("The element is not in the NAO database.")
        return []

def findNAOByElementANDFunctional(element, functional):
    NAOdictionary = openNAODatabase()
    if element in NAOdictionary:
        if functional in NAOdictionary[element]:
            return list(NAOdictionary[element][functional][zetaSelection].keys())
        else:
            print("The functional is not in the NAO database.")
            return []
    else:
        print("The element is not in the NAO database.")
        return []
'''
def findNAOByCutoff(element, threshold, functional = 'gga', zetaSelection = 'dzp' , mode = 'radius_min'):
    '''
    # findNAOByCutoff\n
    @param element: the element of the NAO\n
    @param threshold: the threshold of the radius cutoff\n
    find numerical atomic orbital of one element with specific functional, zeta level in one of the following modes:\n
    1. radius_min: find the NAO with the smallest radius cutoff larger than the threshold\n
    2. radius_max: find the NAO with the largest radius cutoff smaller than the threshold\n
    3. energy_min: find the NAO with the smallest energy cutoff larger than the threshold\n
    4. energy_max: find the NAO with the largest energy cutoff smaller than the threshold\n
    for the following modes, the threshold is a list of two numbers, the first number is the threshold of the energy cutoff, the second number is the threshold of the radius cutoff\n
    5. energy_min_radius_min: find the NAO with the smallest energy cutoff larger than the threshold and the smallest radius cutoff larger than the threshold\n
    6. energy_min_radius_max: find the NAO with the smallest energy cutoff larger than the threshold and the largest radius cutoff smaller than the threshold\n
    7. energy_max_radius_min: find the NAO with the largest energy cutoff smaller than the threshold and the smallest radius cutoff larger than the threshold\n
    8. energy_max_radius_max: find the NAO with the largest energy cutoff smaller than the threshold and the largest radius cutoff smaller than the threshold\n
    @param functional: the functional of the NAO\n
    @param zetaSelection: the zeta notation of the NAO\n
    @param mode: the mode of the search\n
    @return: the NAO file name  
    '''
    NAOdictionary = openNAODatabase()
    if element in NAOdictionary:
        radiusCutoffCache = -1
        NAOfileNameCache = ''
        energyCutoffCache = -1
        if functional in NAOdictionary[element]:
            for fileName in NAOdictionary[element][functional][zetaSelection]:
                if mode == 'radius_min':
                    # for selecting a NAO with the smallest radius cutoff, the radius cutoff must be larger than the threshold
                    # say, add a definitive value to the smallest radius cutoff
                    if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] >= threshold:
                        if NAOfileNameCache == '':
                            radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                            NAOfileNameCache = fileName
                        else:
                            if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] < radiusCutoffCache:
                                radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                                NAOfileNameCache = fileName
                elif mode == 'radius_max':
                    # for selecting a NAO with the largest radius cutoff, the radius cutoff must be smaller than the threshold
                    # say, add a definitive value to the largest radius cutoff
                    if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] <= threshold:
                        if NAOfileNameCache == '':
                            radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                            NAOfileNameCache = fileName
                        else:
                            if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] > radiusCutoffCache:
                                radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                                NAOfileNameCache = fileName

                elif mode == 'energy_min':
                    # similarly.
                    if NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] >= threshold:
                        if NAOfileNameCache == '':
                            energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                            NAOfileNameCache = fileName
                        else:
                            if NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] < energyCutoffCache:
                                energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                                NAOfileNameCache = fileName

                elif mode == 'energy_max':
                    if NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] <= threshold:
                        if NAOfileNameCache == '':
                            energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                            NAOfileNameCache = fileName
                        else:
                            if NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] > energyCutoffCache:
                                energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                                NAOfileNameCache = fileName
 
                elif mode == 'energy_min_radius_min':
                    if NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] >= threshold[0]:
                        if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] >= threshold[1]:
                            if NAOfileNameCache == '':
                                radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                                energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                                NAOfileNameCache = fileName
                            else:
                                if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] < radiusCutoffCache and NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] < energyCutoffCache:
                                    radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                                    energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                                    NAOfileNameCache = fileName
                elif mode == 'energy_min_radius_max':
                    if NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] >= threshold[0]:
                        if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] <= threshold[1]:
                            if NAOfileNameCache == '':
                                radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                                energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                                NAOfileNameCache = fileName
                            else:
                                if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] > radiusCutoffCache and NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] < energyCutoffCache:
                                    radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                                    energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                                    NAOfileNameCache = fileName
                elif mode == 'energy_max_radius_min':
                    if NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] <= threshold[0]:
                        if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] >= threshold[1]:
                            if NAOfileNameCache == '':
                                radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                                energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                                NAOfileNameCache = fileName
                            else:
                                if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] < radiusCutoffCache and NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] > energyCutoffCache:
                                    radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                                    energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                                    NAOfileNameCache = fileName
                elif mode == 'energy_max_radius_max':
                    if NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] <= threshold[0]:
                        if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] <= threshold[1]:
                            if NAOfileNameCache == '':
                                radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                                energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                                NAOfileNameCache = fileName
                            else:
                                if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] > radiusCutoffCache and NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] > energyCutoffCache:
                                    radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                                    energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                                    NAOfileNameCache = fileName
                elif mode == 'energy_min_radius_min':
                    if NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] >= threshold[0]:
                        if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] >= threshold[1]:
                            if NAOfileNameCache == '':
                                radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                                energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                                NAOfileNameCache = fileName
                            else:
                                if NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff'] < radiusCutoffCache and NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff'] < energyCutoffCache:
                                    radiusCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['radiusCutoff']
                                    energyCutoffCache = NAOdictionary[element][functional][zetaSelection][fileName]['energyCutoff']
                                    NAOfileNameCache = fileName
                else:
                    print("Based on present selection rule, numerical orbital is not found in the NAO database.")
            return NAOfileNameCache
        else:
            raise ValueError("The functional is not in the NAO database.")
    else:
        raise ValueError("The element is not in the NAO database.")

# unit test
if __name__ == "__main__":
    print("This module is not intended to be run directly.")
    