import os
import re
import json
from fp2a_keywordsProcessing import keywordsWrite, keywordsRead
def NAODatabaseGenerator(
    path = '.', 
    NAONamePattern = r'[A-Z][a-z]?\w?_[a-z]{3}_([0-9]+)([\.]*)([0-9]*)au_([0-9]*)Ry_([0-9][a-z])+\.orb',
    boolSaveToJsonFile = False,
    jsonFileNameWithPath = './fp2a_NAO.json',
    debugMode = False
    ):
    '''
    This function is used to generate a dictionary of ABACUS NAOs whose name in format of
    [element]_[functional]_[radiusCutoff]au_[energyCutoff]Ry_[valenceElectrons].orb\n
    The dictionary is saved in a json file named fp2a_NAO.json in the current directory by default.\n
    The dictionary is in the format of:\n
    NAOdictionary = {
        element: {
            functional: {
                fileName: {
                    radiusCutoff: radiusCutoff,
                    energyCutoff: energyCutoff,
                    valenceElectrons: valenceElectrons
                    }
                }
            }
        }\n
    @param path: the path of the directory where the NAOs are stored\n
    @param NAONamePattern: the regular expression pattern of the NAO name\n
    @param boolSaveToJsonFile: whether to save the dictionary in a json file\n
    @param jsonFileNameWithPath: the path of the json file\n
    @param debugMode: whether to print the debug information\n
    @return: the dictionary of NAOs
    '''
    NAOdictionary = {}
    fileNames = os.listdir(path)
    for fileName in fileNames:
        match = re.search(NAONamePattern, fileName)
        if match:
            words = re.split(pattern = "_|.orb|au|Ry", string = fileName)
            element = keywordsRead(words[0])
            functional = keywordsRead(words[1])
            radiusCutoff = keywordsRead(words[2].replace('au', ''))
            energyCutoff = keywordsRead(words[4].replace('Ry', ''))
            valenceElectrons = keywordsRead(words[6])
            if debugMode:
                print("fileName: ", keywordsWrite(fileName))
                print("element: ", keywordsWrite(element))
                print("functional: ", keywordsWrite(functional))
                print("radiusCutoff: ", keywordsWrite(radiusCutoff))
                print("energyCutoff: ", keywordsWrite(energyCutoff))
                print("valenceElectrons: ", keywordsWrite(valenceElectrons))
            if element not in NAOdictionary:
                NAOdictionary[element] = {}
            if functional not in NAOdictionary[element]:
                NAOdictionary[element][functional] = {}
            NAOdictionary[element][functional][fileName] = {
                "radiusCutoff": radiusCutoff,
                "energyCutoff": energyCutoff,
                "valenceElectrons": valenceElectrons
                }
    if boolSaveToJsonFile:
        with open(jsonFileNameWithPath, 'w') as fp:
            json.dump(NAOdictionary, fp, indent = 4)
    return NAOdictionary

def openNAODatabase():
    if os.path.isfile('./fp2a_NAO.json'):
        with open('./fp2a_NAO.json', 'r') as fp:
            NAOdictionary = json.load(fp)
        return NAOdictionary
    else:
        return NAODatabaseGenerator()

def findNAOByElement(element):
    NAOdictionary = openNAODatabase()
    if element in NAOdictionary:
        returnList = []
        for functional in NAOdictionary[element]:
            returnList += list(NAOdictionary[element][functional].keys())
        return returnList
    else:
        print("The element is not in the NAO database.")
        return []

def findNAOByElementANDFunctional(element, functional):
    NAOdictionary = openNAODatabase()
    if element in NAOdictionary:
        if functional in NAOdictionary[element]:
            return list(NAOdictionary[element][functional].keys())
        else:
            print("The functional is not in the NAO database.")
            return []
    else:
        print("The element is not in the NAO database.")
        return []

if __name__ == "__main__":
    pass
    