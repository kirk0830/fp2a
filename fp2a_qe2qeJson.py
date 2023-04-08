from re import T, split as resplit
import json
import numpy as np
# pre-set data
nonStandardSectionTitles = [
    "ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS",
    "ADDITIONAL_K_POINTS", "CELL_PARAMETERS", "CONSTRAINTS",
    "OCCUPATIONS", "ATOMIC_VELOCITIES", "ATOMIC_FORCES",
    "SOLVENTS", "HUBBARD"
    ]
nSSTwithSectionParameter = [
    "ATOMIC_POSITIONS", "K_POINTS", 
    "ADDITIONAL_K_POINTS", "CELL_PARAMETERS",
    "HUBBARD", "SOLVENTS"
]
# =================== functions ===================
def keywordsParse(keywordDefinition):
    # return a dict of keyword and value
    equalIndex = keywordDefinition.find('=')
    if equalIndex != -1:
        # typical keyword definition line
        # e.g., calculation = 'scf'
        words = keywordDefinition.replace(' ', '').replace('\'','').split('=')
        return words
    else:
        newDefinition = keywordDefinition
        while (newDefinition.find("  ") != -1):
            newDefinition = newDefinition.replace("  ", " ")
        spaceIndex = newDefinition.find(" ")
        return [newDefinition[0:spaceIndex], newDefinition[spaceIndex+1::]]
def keywordPostProcessing(keyword):
    
    keywordLowerCase = keyword.lower()
    if keywordLowerCase == ".true.":
        return True
    elif keywordLowerCase == ".false.":
        return False
    else:
        dotIndex = keyword.find('.')
        if dotIndex != -1:
            try:
                return float(keyword)
            except ValueError:
                scientificNotation = resplit('([eEdD])', keyword)
                try:
                    return float(scientificNotation[0]) * 10 ** int(scientificNotation[2])
                except ValueError:
                    return keyword
        else:
            try:
                return int(keyword)
            except ValueError:
                return keyword
def isNonStandardSectionTitle(line):
    words = line.strip().split()
    if nonStandardSectionTitles.count(words[0]) > 0:
        return True
    else:
        return False
def isNSSTwithSectionParameter(line):
    words = line.strip().split()
    if nSSTwithSectionParameter.count(words[0]) > 0:
        return True
    else:
        return False

def qe2qeJson(fileName, boolVersion7_1 = True, parseMode = 'single', boolDebugMode = False):

    """
    # Quantum ESPRESSO input script parser\n
    @param fileName: the file name of the input script\n
    @param boolVersion7_1: whether the input script is for QE version >= 7.1\n
    @param parseMode: 'single' or 'examples', 'single' means the input script is a single input script, 'examples' means the input script is a run_example file which contains multiple input scripts\n
    @param boolDebugMode: whether to print debug information\n
    @return: a list of dictionary-organized quantum espresso input script(s)\n
    @note: the input script can be run_example file which contains not-only one input script,
    but certainly the most common kind of input scripts is supported. (laugh)
    """
    returnList = []
    quantumESPRESSOInput = {}

    boolReadInputScript = False
    boolReadStandardSection = False
    boolReadNonstandardSection = False
    boolReadAtomicSpeciesSection = False
    boolReadAtomicPositionsSection = False
    boolReadCellParametersSection = False
    boolReadKpointsSection = False
    boolReadKpath = False
    boolReadHubbardSection = False # supported only for qe versions >= 7.1
    boolReadSolventsSection = False # supported only for qe versions >= 7.1
    intKpointSectionLineRead = 0
    intKpointSectionLineToRead = -1
    sectionTitle = "__NONE__"

    def hubbardSectionStartswith(line):
        return line.startswith("U ") or line.startswith("J ") or line.startswith("V ") or line.startswith("B ")

    def checkStatus():
        print("CHECK STATUS:")
        print("boolReadInputScript: "+str(boolReadInputScript))
        print("boolReadStandardSection: "+str(boolReadStandardSection))
        print("boolReadNonstandardSection: "+str(boolReadNonstandardSection))
        print("boolReadAtomicSpeciesSection: "+str(boolReadAtomicSpeciesSection))
        print("boolReadAtomicPositionsSection: "+str(boolReadAtomicPositionsSection))
        print("boolReadCellParametersSection: "+str(boolReadCellParametersSection))
        print("boolReadKpointsSection: "+str(boolReadKpointsSection))
        print("boolReadKpath: "+str(boolReadKpath))
        print("boolVersion7_1: "+str(boolVersion7_1))
        print("boolReadHubbardSection: "+str(boolReadHubbardSection))
        print("boolReadSolventsSection: "+str(boolReadSolventsSection))
        print("intKpointSectionLineRead: "+str(intKpointSectionLineRead))
        print("intKpointSectionLineToRead: "+str(intKpointSectionLineToRead))
        print("sectionTitle: "+sectionTitle)
    
    def resetStatus():
        boolReadInputScript = False
        boolReadStandardSection = False
        boolReadNonstandardSection = False
        boolReadAtomicSpeciesSection = False
        boolReadAtomicPositionsSection = False
        boolReadCellParametersSection = False
        boolReadKpointsSection = False
        boolReadKpath = False
        boolReadHubbardSection = False # supported only for qe versions >= 7.1
        boolReadSolventsSection = False # supported only for qe versions >= 7.1
        intKpointSectionLineRead = 0
        intKpointSectionLineToRead = -1
        sectionTitle = "__NONE__"

    def createHubbardSection():

        if not boolVersion7_1:
            try:
                if quantumESPRESSOInput["system"]["lda_plus_u"]:
                    quantumESPRESSOInput["HUBBARD"] = {
                        "on-site": {},
                        "hopping": {}
                    }
                    # manifold information?
                    for i, element in enumerate(quantumESPRESSOInput["ATOMIC_SPECIES"]["elements"]):
                        quantumESPRESSOInput["HUBBARD"]["on-site"][element] = {
                            "unknown_manifold_from_old_qe_versions":{
                                "U": quantumESPRESSOInput["system"]["Hubbard_U({0})".format(str(i+1))]
                            }
                        }
            except KeyError:
                # no lda_plus_u
                pass
        else:
            pass

    def createCellParameterSection():
        if "CELL_PARAMETERS" not in quantumESPRESSOInput.keys():
            if quantumESPRESSOInput["system"]["ibrav"] == 0:
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                    }
            elif quantumESPRESSOInput["system"]["ibrav"] == 1:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [[a, 0, 0], [0, a, 0], [0, 0, a]]
                    }
            elif quantumESPRESSOInput["system"]["ibrav"] == 2:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [[-a/2, 0, a/2], [0, a/2, a/2], [-a/2, a/2, 0]]
                    }
            elif quantumESPRESSOInput["system"]["ibrav"] == 3:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [[a/2, a/2, a/2], [-a/2, a/2, a/2], [-a/2, -a/2, a/2]]
                    }
            elif quantumESPRESSOInput["system"]["ibrav"] == -3:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [[-a/2, a/2, a/2], [a/2, -a/2, a/2], [a/2, a/2, -a/2]]
                    }
            elif quantumESPRESSOInput["system"]["ibrav"] == 4:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [[a, 0, 0], [-a/2, np.sqrt(3)/2*a, 0], [0, 0, c]]
                    }
            elif quantumESPRESSOInput["system"]["ibrav"] == 5:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                cosgamma = quantumESPRESSOInput["system"]["celldm(4)"]
                tx = np.sqrt((1-cosgamma)/2)
                ty = np.sqrt((1-cosgamma)/6)
                tz = np.sqrt((1+2*cosgamma)/3)
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [[a*tx, -a*ty, a*tz], [0, 2*a*ty, a*tz], [-a*tx, -a*ty, a*tz] ]
                    }
            elif quantumESPRESSOInput["system"]["ibrav"] == -5:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                cosgamma = quantumESPRESSOInput["system"]["celldm(4)"]
                tx = np.sqrt((1-cosgamma)/2)
                ty = np.sqrt((1-cosgamma)/6)
                tz = np.sqrt((1+2*cosgamma)/3)
                u = tz - 2*np.sqrt(2)*ty
                v = tz + np.sqrt(2)*ty
                a_prime = a/np.sqrt(3)

                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [
                        [a_prime*u, a_prime*v, a_prime*v],
                        [a_prime*v, a_prime*u, a_prime*v],
                        [a_prime*v, a_prime*v, a_prime*u],
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == 6:
                a = quantumESPRESSOInput["system"]["celldm(1)"] 
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [[a, 0, 0], [0, a, 0], [0, 0, c]
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == 7:
                a = quantumESPRESSOInput["system"]["celldm(1)"] 
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [[a/2, -a/2, c/2], [a/2, a/2, c/2], [-a/2, -a/2, c/2],
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == 8:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                b = quantumESPRESSOInput["system"]["celldm(2)"]*a
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [[a, 0, 0], [0, b, 0], [0, 0, c]
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == 9:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                b = quantumESPRESSOInput["system"]["celldm(2)"]*a
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [
                        [a/2, b/2, 0], [-a/2, b/2, 0], [0, 0, c]
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == -9:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                b = quantumESPRESSOInput["system"]["celldm(2)"]*a
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [
                        [a/2, -b/2, 0], [a/2, b/2, 0], [0, 0, c]
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == 91:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                b = quantumESPRESSOInput["system"]["celldm(2)"]*a
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [
                        [a, 0, 0], [0, b/2, -c/2], [0, b/2, c/2]
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == 10:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                b = quantumESPRESSOInput["system"]["celldm(2)"]*a
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [
                        [a/2, 0, c/2], [a/2, b/2, 0], [0, b/2, c/2]
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == 11:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                b = quantumESPRESSOInput["system"]["celldm(2)"]*a
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [
                        [a/2, b/2, c/2],
                        [-a/2, b/2, c/2],
                        [-a/2, -b/2, c/2]
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == 12:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                b = quantumESPRESSOInput["system"]["celldm(2)"]*a
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                cos_gamma = quantumESPRESSOInput["system"]["celldm(4)"]
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [
                        [a, 0, 0],
                        [b*cos_gamma, b*np.sqrt(1-cos_gamma**2), 0],
                        [0, 0, c]
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == -12:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                b = quantumESPRESSOInput["system"]["celldm(2)"]*a
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                cos_ac = quantumESPRESSOInput["system"]["celldm(5)"]
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [
                        [a, 0, 0],
                        [0, b, 0],
                        [c*cos_ac, 0, c*np.sqrt(1-cos_ac**2)]
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == 13:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                b = quantumESPRESSOInput["system"]["celldm(2)"]*a
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                cos_gamma = quantumESPRESSOInput["system"]["celldm(4)"]
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [
                        [a/2, 0, -c/2],
                        [b*cos_gamma, b*np.sqrt(1-cos_gamma**2), 0],
                        [a/2, 0, c/2]
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == -13:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                b = quantumESPRESSOInput["system"]["celldm(2)"]*a
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                cos_beta = quantumESPRESSOInput["system"]["celldm(5)"]
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [
                        [a/2, b/2, 0],
                        [-a/2, b/2, 0],
                        [c*cos_beta, 0, c*np.sqrt(1-cos_beta**2)]
                    ]
                }
            elif quantumESPRESSOInput["system"]["ibrav"] == 14:
                a = quantumESPRESSOInput["system"]["celldm(1)"]
                b = quantumESPRESSOInput["system"]["celldm(2)"]*a
                c = quantumESPRESSOInput["system"]["celldm(3)"]*a
                cos_alpha = quantumESPRESSOInput["system"]["celldm(4)"]
                cos_beta = quantumESPRESSOInput["system"]["celldm(5)"]
                cos_gamma = quantumESPRESSOInput["system"]["celldm(6)"]
                sin_gamma = np.sqrt(1-cos_gamma**2)
                v1 = [a, 0, 0]
                v2 = [b*cos_gamma, b*sin_gamma, 0]
                v3 = [c*cos_beta, 
                    c*(cos_alpha-cos_beta*cos_gamma)/sin_gamma,
                    c*np.sqrt(1+2*cos_alpha*cos_beta*cos_gamma-cos_alpha**2-cos_beta**2-cos_gamma**2)/sin_gamma]
                quantumESPRESSOInput["CELL_PARAMETERS"] = {
                    "unit": "angstrom",
                    "cell": [
                        v1, v2, v3
                    ]
                }
                                  

    with open(file = fileName, mode = 'r') as f:
        
        lines = f.readlines()

        for idx, line in enumerate(lines):
            cleanLine = line.replace('\n','').replace('\'','').strip()
            if len(cleanLine) == 0: continue
            if boolReadInputScript and boolDebugMode:
                print("PARSE_LINE_STARTS->|"+cleanLine+"|<-END")
                pass
            # =================== read section titles ===================
            if cleanLine.startswith('&'):
                if cleanLine.startswith('&control') or cleanLine.startswith('&CONTROL'):
                    boolReadInputScript = True
                    if boolVersion7_1:
                        quantumESPRESSOInput['__version__'] = '>=7.1'
                    else:
                        quantumESPRESSOInput['__version__'] = '<7.1'
                else:
                    if not boolReadInputScript:
                        print("current line read: "+cleanLine)
                        print("present input script is not organized in common sequence? not control section be the first? check it")
                # qe real input script read from here.
                boolReadStandardSection = True
                sectionTitle = cleanLine[1::].lower()
                quantumESPRESSOInput[sectionTitle] = {}
            elif isNonStandardSectionTitle(cleanLine):
                boolReadNonstandardSection = True
                sectionTitle = cleanLine
                unit = 'absent_here_so_it_should_be_default'
                if isNSSTwithSectionParameter(cleanLine):
                    words = cleanLine.split()
                    if len(words) > 1:
                        unit = words[1].replace('(','').replace(')','').replace('{','').replace('}','')
                        sectionTitle = cleanLine.split()[0]
                if sectionTitle == "K_POINTS":
                    if unit.lower() == 'gamma':
                        quantumESPRESSOInput["K_POINTS"] = {"mode": "mk", "unit": unit,"grid":[1, 1, 1], "shift":[0, 0, 0]}
                        boolReadKpath = False
                        boolReadKpointsSection = False
                        boolReadNonstandardSection = False
                    else:
                        boolReadKpointsSection = True
                        theNextLineSplit = lines[idx+1].strip().split()
                        if len(theNextLineSplit) < 6:
                            boolReadKpath = True
                            quantumESPRESSOInput["K_POINTS"] = {"mode": "kpath", "unit": unit, "kpts":[], "w_kpts":[]}
                            intKpointSectionLineToRead = int(theNextLineSplit[0])
                        else:
                            quantumESPRESSOInput["K_POINTS"] = {"mode": "mk", "unit": unit,"grid":[], "shift":[]}
                            boolReadKpath = False
                elif sectionTitle == "ATOMIC_SPECIES":
                    boolReadAtomicSpeciesSection = True
                    quantumESPRESSOInput["ATOMIC_SPECIES"] = {"elements":[], "masses":[], "pseudopotentials":[]}
                elif sectionTitle == "ATOMIC_POSITIONS":
                    boolReadAtomicPositionsSection = True
                    quantumESPRESSOInput["ATOMIC_POSITIONS"] = {"unit": unit, "elements":[], "coordinates":[], "constraints":[]}
                elif sectionTitle == "CELL_PARAMETERS":
                    boolReadCellParametersSection = True
                    quantumESPRESSOInput["CELL_PARAMETERS"] = {"unit": unit,"cell":[]}

                elif sectionTitle == "HUBBARD" and boolVersion7_1:
                    boolReadHubbardSection = True
                    theNextLineSplit = lines[idx+1].strip().split()
                    quantumESPRESSOInput["HUBBARD"] = {"on-site": {}, "hopping": {}}
                    # the format of Hubbard section is suggested to be organized in the way above.
                    # generally, there will be two kinds of dataline in this section,
                    # the first like: U Fe-3d 5.0, length = 3
                    # the second like: J Fe-3d O-2p 1 15 0.0, length = 6
                    # and the data will save in this way:
                    # {
                    # 'on-site': {'Fe': 
                    #               {
                    #               '4s': {'U': 2.0, 'J': 0.5, 'B': 0.0}, 
                    #               '3d': {'U': 2.0, 'J': 0.5, 'B': 0.0}
                    #               }
                    #             'O':...
                    #             }, 
                    # 'hopping': {'1-4s-19-2p': ['Fe', 'O', 0.75]}
                    # }
                elif sectionTitle == "SOLVENTS" and boolVersion7_1:
                    boolReadSolventsSection = True
                    quantumESPRESSOInput["SOLVENTS"] = {"unit": unit, "labels":[], "densities":[], "molecules":[]}
                else:
                    quantumESPRESSOInput[sectionTitle] = {}
            # =================== non-title lines =======================
            else:
                if not boolReadInputScript:
                    # linux command line
                    continue
                else:
            # ============= keywords of standard section ================
                    if cleanLine.startswith('!'): continue # comment line
                    if boolReadStandardSection:
                        if cleanLine.startswith('/'):
                            boolReadStandardSection = False
                        else:
                            # enter standard section and parse keyword definitions
                            commaIndex = cleanLine.find(',')
                            if commaIndex != -1:
                                # meaning there are keywords definition seperated by comma
                                # this happens usually in standard section
                                keywordsDefinition = cleanLine.split(',')
                                keywordsDefinition = [
                                    keywordDefinition.strip() for keywordDefinition in keywordsDefinition if len(keywordDefinition.strip()) > 0
                                    ]
                                for keywordDefinition in keywordsDefinition:
                                    keywordsKeyValuePair = keywordsParse(keywordDefinition)
                                    quantumESPRESSOInput[sectionTitle][keywordsKeyValuePair[0]] = keywordPostProcessing(keywordsKeyValuePair[1])
                            else:
                                # quite common word definition line
                                keywordsKeyValuePair = keywordsParse(cleanLine)
                                quantumESPRESSOInput[sectionTitle][keywordsKeyValuePair[0]] = keywordPostProcessing(keywordsKeyValuePair[1])
            # ============= keywords of non-standard section ============
            # not standard section, e.g., ATOMIC_SPECIES
            # but addtionally for KPOINTS, keywords should not be read as other sections
            # for case reading kpoints line
                    else:
                        if boolReadNonstandardSection:
                            words = cleanLine.split()
                            if boolReadKpointsSection:
                                if boolReadKpath:
                                    if len(words) < 2:
                                        #print("It should be the number of kpts read: "+cleanLine)
                                        try:
                                            intKpointSectionLineToRead = int(cleanLine)
                                            #print("The number of kpts read is correct")
                                        except ValueError:
                                            checkStatus()
                                            exit(1)
                                    else:
                                        intKpointSectionLineRead += 1
                                        if intKpointSectionLineRead < intKpointSectionLineToRead:
                                            quantumESPRESSOInput["K_POINTS"]["kpts"].append([
                                                keywordPostProcessing(words[0]), keywordPostProcessing(words[1]), keywordPostProcessing(words[2])])
                                            quantumESPRESSOInput["K_POINTS"]["w_kpts"].append(keywordPostProcessing(words[3]))
                                        else:
                                            boolReadKpath = False
                                            boolReadKpointsSection = False
                                            boolReadNonstandardSection = False
                                else:
                                    quantumESPRESSOInput["K_POINTS"]["grid"] = [
                                        keywordPostProcessing(words[0]), 
                                        keywordPostProcessing(words[1]), 
                                        keywordPostProcessing(words[2])
                                        ]
                                    quantumESPRESSOInput["K_POINTS"]["shift"] = [
                                        keywordPostProcessing(words[3]), 
                                        keywordPostProcessing(words[4]), 
                                        keywordPostProcessing(words[5])
                                        ]
                                    boolReadKpointsSection = False
                                    boolReadNonstandardSection = False
                            elif boolReadAtomicSpeciesSection:
                                quantumESPRESSOInput["ATOMIC_SPECIES"]["elements"].append(words[0])
                                quantumESPRESSOInput["ATOMIC_SPECIES"]["masses"].append(keywordPostProcessing(words[1]))
                                quantumESPRESSOInput["ATOMIC_SPECIES"]["pseudopotentials"].append(words[2])
                                if len(quantumESPRESSOInput["ATOMIC_SPECIES"]["elements"]) == quantumESPRESSOInput["system"]["ntyp"]:
                                    boolReadAtomicSpeciesSection = False
                                    boolReadNonstandardSection = False
                            elif boolReadAtomicPositionsSection:
                                quantumESPRESSOInput["ATOMIC_POSITIONS"]["elements"].append(words[0])
                                quantumESPRESSOInput["ATOMIC_POSITIONS"]["coordinates"].append([
                                    keywordPostProcessing(words[1]), keywordPostProcessing(words[2]), keywordPostProcessing(words[3])])
                                if len(words) == 7:
                                    quantumESPRESSOInput["ATOMIC_POSITIONS"]["constraints"].append([
                                        keywordPostProcessing(words[4]), keywordPostProcessing(words[5]), keywordPostProcessing(words[6])])
                                else:
                                    quantumESPRESSOInput["ATOMIC_POSITIONS"]["constraints"].append([1, 1, 1])
                                if len(quantumESPRESSOInput["ATOMIC_POSITIONS"]["elements"]) == quantumESPRESSOInput["system"]["nat"]:
                                    boolReadAtomicPositionsSection = False
                                    boolReadNonstandardSection = False
                            elif boolReadCellParametersSection:
                                quantumESPRESSOInput["CELL_PARAMETERS"]["cell"].append([
                                    keywordPostProcessing(words[0]), 
                                    keywordPostProcessing(words[1]), 
                                    keywordPostProcessing(words[2])
                                    ])
                                if len(quantumESPRESSOInput["CELL_PARAMETERS"]["cell"]) == 3:
                                    boolReadCellParametersSection = False
                                    boolReadNonstandardSection = False
                            elif boolReadHubbardSection and boolVersion7_1:
                                # for nonstandard section, there is no signal where it needs, but one needs to check manually
                                # {
                                # 'on-site': {'Fe': 
                                #               {
                                #               '4s': {'U': 2.0, 'J': 0.5, 'B': 0.0}, 
                                #               '3d': {'U': 2.0, 'J': 0.5, 'B': 0.0}
                                #               }
                                #             'O':...
                                #             }, 
                                # 'hopping': {'1-4s-19-2p': ['Fe', 'O', 0.75]}
                                # }
                                if hubbardSectionStartswith(cleanLine):
                                    if len(words) == 3:
                                        # U Fe-3d 0.5
                                        specificInformation = words[1].split("-")
                                        element = specificInformation[0]
                                        manifold = specificInformation[1]
                                        if element not in quantumESPRESSOInput["HUBBARD"]["on-site"]:
                                            quantumESPRESSOInput["HUBBARD"]["on-site"][element] = {}
                                        if manifold not in quantumESPRESSOInput["HUBBARD"]["on-site"][element]:
                                            quantumESPRESSOInput["HUBBARD"]["on-site"][element][manifold] = {}
                                        quantumESPRESSOInput["HUBBARD"]["on-site"][element][manifold][words[0]] = keywordPostProcessing(words[2])
                                    elif len(words) == 6:
                                        # V Fe-3d O-2p 1 19 0.75
                                        specificInformation1 = words[1].split("-")
                                        element1 = specificInformation1[0]
                                        manifold1 = specificInformation1[1]
                                        atomIndex1 = words[3]
                                        specificInformation2 = words[2].split("-")
                                        element2 = specificInformation2[0]
                                        manifold2 = specificInformation2[1]
                                        atomIndex2 = words[4]
                                        quantumESPRESSOInput["HUBBARD"]["hopping"]["-".join([atomIndex1, manifold1, atomIndex2, manifold2])] = [
                                            element1, element2, keywordPostProcessing(words[5])]
                                else:
                                    # it is the way to leave HUBBARD the non standard section
                                    if cleanLine == "EOF" and parseMode == 'examples':
                                        if len(quantumESPRESSOInput.keys()) > 0:
                                            createHubbardSection()
                                            createCellParameterSection()
                                            resetStatus()
                                            returnList.append(quantumESPRESSOInput)
                                            quantumESPRESSOInput = {}
                                    # there is a bug of this program: for single mode if HUBBARD section is not the last section
                                    elif parseMode == 'single':
                                        createHubbardSection()
                                        resetStatus()
                            elif boolReadSolventsSection and boolVersion7_1:
                                quantumESPRESSOInput["SOLVENTS"]["labels"].append(words[0])
                                quantumESPRESSOInput["SOLVENTS"]["densities"].append(keywordPostProcessing(words[1]))
                                quantumESPRESSOInput["SOLVENTS"]["molecules"].append(words[2])
                                if len(quantumESPRESSOInput["SOLVENTS"]["labels"]) == quantumESPRESSOInput["rism"]["nsolv"]:
                                    boolReadSolventsSection = False
                                    boolReadNonstandardSection = False
                            else:
                                keywordsKeyValuePair = keywordsParse(cleanLine)
                                quantumESPRESSOInput[sectionTitle][keywordsKeyValuePair[0]] = keywordsKeyValuePair[1]
                                boolReadNonstandardSection = False
                        else:
                            if cleanLine == "EOF":
                                if len(quantumESPRESSOInput.keys()) > 0:
                                    createHubbardSection()
                                    createCellParameterSection()
                                    resetStatus()
                                    returnList.append(quantumESPRESSOInput)
                                    quantumESPRESSOInput = {}

            if parseMode == 'single' and idx == len(lines) - 1:
                createHubbardSection()
                createCellParameterSection()
                resetStatus()
                if len(quantumESPRESSOInput.keys()) > 0:
                    returnList.append(quantumESPRESSOInput)
                    quantumESPRESSOInput = {}

    return returnList

if __name__ == "__main__":


    inputScriptFullFileName = "EXX_example\\run_example"
    jsonFileNamePrefix = inputScriptFullFileName.replace('\\', '_').replace('.', '_')
    DictionariesList = qe2qeJson(inputScriptFullFileName, boolVersion7_1 = True, parseMode='examples')
    for i, dictionary in enumerate(DictionariesList):
        #print("Dictionary "+str(i)+": "+str(dictionary))
        with open(jsonFileNamePrefix +'_'+ str(i) + '.json', 'w') as outfile:
            json.dump(dictionary, outfile, indent=4)

    """
    current_dir = osgetcwd()
    print("Validity check for input script parsing")
    print("Current directory: "+current_dir)
    directories = next(oswalk(current_dir))[1]
    for directory in directories:
        print("Test runs in directory: "+directory)
        files = [f for f in oslistdir(directory) if not f.startswith('.')]
        for file in files:
            if file.startswith("run_example"):
                print("Test run for file: "+file)
                qe2qeJson(f"{directory}\\run_example")
    """
