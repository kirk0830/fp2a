import json
# only for rotating dictionary of specific format
def jsonRotate(rotateSequence, sourceDictionaryFileName = "depuis_le_abacus.json"):

    rotatedDictionaryFileName = sourceDictionaryFileName[:-5] + "_rotated-"+rotateSequence+".json"
    with open(sourceDictionaryFileName) as sourceDictionaryFile:
        sourceDictionary = json.load(sourceDictionaryFile)["corresponding"]
        rotatedDictionary = {}
        for key in sourceDictionary:
            if rotateSequence == "qac":
                rotatedDictionary[sourceDictionary[key][0]] = [key, sourceDictionary[key][1], sourceDictionary[key][2]]
            elif rotateSequence == "cqa":
                rotatedDictionary[sourceDictionary[key][1]] = [sourceDictionary[key][0], key, sourceDictionary[key][2]]
            elif rotateSequence == "acq":
                rotatedDictionary[key] = [sourceDictionary[key][1], sourceDictionary[key][0], sourceDictionary[key][2]]
            elif rotateSequence == "aqc":
                rotatedDictionary[key] = [sourceDictionary[key][0], sourceDictionary[key][1], sourceDictionary[key][2]]
            elif rotateSequence == "caq":
                rotatedDictionary[sourceDictionary[key][1]] = [key, sourceDictionary[key][0], sourceDictionary[key][2]]
            elif rotateSequence == "qca":
                rotatedDictionary[sourceDictionary[key][0]] = [sourceDictionary[key][1], key, sourceDictionary[key][2]]
        rotatedDictionary.pop("NaN")
        
    with open(rotatedDictionaryFileName, 'w') as rotatedDictionaryFile:
        json.dump(rotatedDictionary, rotatedDictionaryFile, indent=4)
    return rotatedDictionaryFileName

if __name__ == '__main__':
    rotateSequence = "qac" # means qe keyword is the key and abacus, cp2k keywords are the values
    jsonRotate(rotateSequence)