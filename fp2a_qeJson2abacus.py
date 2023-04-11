import json
import os
from fp2a_keywordsProcessing import keywordsWrite
from fp2a_qe2qeJson import isNonStandardSectionTitle
from fp2a_numericalOrbitalManager import findNAOByCutoff
from fp2a_jsonRotation import jsonRotate

def discrepantKeywordsConversion(keyword, value):
    # This function can only complete single line conversion
    # for more complicated version, program in generateABACUSFiles function :(
    if keyword == 'calculation':
        optionListDirect = ['scf', 'relax', 'nscf', 'md']
        if optionListDirect.count(value) > 0:
            return value
        elif value == 'vc-relax':
            return 'cell-relax'        
        else:
            print('|-Conversion warning: ' + value + ' is not a valid value for keyword ' + keyword + ', using default value')
            return 'scf'
    if keyword == 'ibrav':
        ibravReturnDict = {
            "0": "none", "1": "cubic", "2": "fcc",  "3": "bcc", "4": "hexagonal", "5": "triangoal", "6": "st",
            "7": "bct", "8": "so", "9": "baco", "10": "fco", "11": "bco", "12": "sm", "13": "bacm", "14": "triclinic"
        }
        if str(value) in ibravReturnDict:
            return ibravReturnDict[str(value)]
        else:
            raise ValueError('*ERROR: ibrav = ' + str(value) + ' is not supported in ABACUS')
    else:
        return value

def discrepantModuleConversion(qeDict, module = 'smearing'):


    if module == 'smearing':
        print('|-Runtime information: keyword \'occupations\' is needed to specially treated:')
        if qeDict["system"]["occupations"] == 'smearing':
            print('|-Runtime information: keyword \'occupations\' is set to \'smearing\'')
            if 'smearing' in qeDict["system"]:
                print('|-Runtime information: keyword \'smearing\' is set to \'' 
                    + qeDict["system"]["smearing"] + '\'')
                optionListGaussian = ['gaussian', 'gauss']
                optionListFermiDirac = ['fermi-dirac', 'f-d', 'fd']
                optionListMethfesselPaxton = ['methfessel-paxton', 'm-p', 'mp']
                if optionListGaussian.count(qeDict["system"]["smearing"]) > 0:
                    return ['smearing_method'], ['gaussian']
                elif optionListFermiDirac.count(qeDict["system"]["smearing"]) > 0:
                    return ['smearing_method'], ['fermi-dirac']
                elif optionListMethfesselPaxton.count(qeDict["system"]["smearing"]) > 0:
                    return ['smearing_method'], ['methfessel-paxton']
                else:
                    print('|-Runtime information: keyword \'smearing\' is not supported by ABACUS, set to default \'gaussian\'')
                    return ['smearing_method'], ['gaussian']
            else:
                print('|-Runtime information: keyword \'smearing\' is not set, default is \'gaussian\'')
                return ['smearing_method'], ['gaussian']
        else:
            print('|-Runtime information: keyword \'occupations\' is not supported by ABACUS, set to default \'fixed\'')
            return ['smearing_method'], ['fixed']
        
    elif module == 'vdw':
        optionListD2 = ['grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d']
        optionListD3 = ['grimme-d3', 'Grimme-D3', 'DFT-D3', 'dft-d3']
        if optionListD2.count(qeDict["system"]["vdw_corr"]) > 0:
            return ['vdw_method'], ['d2']
        elif optionListD3.count(qeDict["system"]["vdw_corr"]) > 0:
            if 'dftd3_version' in qeDict["system"]:
                if qeDict["system"]["dftd3_version"] == 2:
                    raise Exception('*ERROR: parameter conflict! you set vdw_corr to dft-d3 but dftd3_version to 2, which means dft-d2')
                elif qeDict["system"]["dftd3_version"] == 3:
                    return ['vdw_method'], ['d3_0']
                elif qeDict["system"]["dftd3_version"] == 4:
                    return ['vdw_method'], ['d3_bj']
                elif qeDict["system"]["dftd3_version"] == 5:
                    print('|-Conversion warning: DFTD3 version 5 (Grimme-D3M (zero damping)) is not implemented in abacus, will use version 3 instead')
                    return ['vdw_method'], ['d3_0']
                elif qeDict["system"]["dftd3_version"] == 6:
                    print('|-Conversion warning: DFTD3 version 6 (Grimme-D3M (BJ damping)) is not implemented in abacus, will use version 4 instead')
                    return ['vdw_method'], ['d3_bj']
                else:
                    raise Exception('*ERROR: DFTD3 version is not valid')
            else:
                return ['vdw_method'], ['d3_0']
        else:
            return ['vdw_method'], ['none']
    #elif module == 'SOLVANTS':
    #    pass
    else:
        return [''], ['']
def hubbardManifoldToNumber(manifold):
    if manifold == 'unknown_manifold_from_old_qe_versions':
        print('|-Hubbard warning: For you are using Quantum ESPRESSO < 7.1, \n'
             +'|                  no information for determining the angular momentum of the element, \n'
             +'|                  therefore ASSUMING DFT+U DISABLED')
        return '-1'
    elif manifold[-1] == 's':
        return '0'
    elif manifold[-1] == 'p':
        return '1'
    elif manifold[-1] == 'd':
        return '2'
    elif manifold[-1] == 'f':
        return '3'
    elif manifold[-1] == 'g':
        return '4'
    elif manifold[-1] == 'h':
        return '5'
    else:
        # are you serious you are using such a high angular momentum?
        return '-1'

def fromElementFindLinesInXYZ(element, elementListInXYZ):

    lines = []
    for i in range(len(elementListInXYZ)):
        if elementListInXYZ[i] == element:
            lines.append(i)
    return lines

def generateABACUSFiles(
    mode = 'on-the-fly',
    qeInputScriptJson = None,
    qeInputScriptJsonFile = 'qeInputScript.json', 
    keywordConversionJson = None,
    keywordConversionJsonFile = 'keywordConversion.json', 
    ABACUSInputFile = 'INPUT', 
    ABACUSStructureFile = 'STRU', 
    ABACUSKpointsFile = 'KPT',
    qeVersion_7_1 = True,
    overwriteKeywords = {},
    additionalKeywords = {"numerical_orbitals": {"use_nao": False}},
    boolUpdateInformation = False
    ):

    """
    # generate ABACUS input files from QE input script in json format\n
    @param mode: the mode controlling how to read dictionary data, available choices are on-the-fly and file\n
    @param qeInputScriptJsonFile: QE input script in json format generated by qe2qeJson.py\n
    @param keywordConversionJsonFile: json file containing the conversion table of QE keywords to ABACUS keywords,
    generated by jsonRotation.py\n
    @param ABACUSInputFile: ABACUS input file name, default is INPUT\n
    @param ABACUSStructureFile: ABACUS structure file name, default is STRU\n
    @param ABACUSKpointsFile: ABACUS kpoints file name, default is KPT\n
    @param qeVersion_7_1: whether the QE version is >= 7.1, default is True\n
    @param overwriteKeywords: keywords to be overwritten, default is empty\n
    @param additionalKeywords: additional keywords to be added to the ABACUS input file, default contains a
    dictionary named numerical orbitals\n
    """
    qeNonstandardSections = []

    if mode == 'on-the-fly':
        if qeInputScriptJson == None:
            qeInputScript = {}
            raise Exception('*ERROR: qeInputScriptJson is not provided')
        if keywordConversionJson == None:
            keywordConversion = {}
            raise Exception('*ERROR: keywordConversionJson is not provided')
        else:
            qeInputScript = qeInputScriptJson
            keywordConversion = keywordConversionJson

    elif mode == 'file':
        if os.path.isfile(qeInputScriptJsonFile) == False:
            raise Exception('*ERROR: qeInputScriptJsonFile does not exist')
        else:
            with open(qeInputScriptJsonFile, 'r') as f:
                qeInputScript = json.load(f)
        if os.path.isfile(keywordConversionJsonFile) == False:
            print('*WARNING: keywordConversionJsonFile does not exist, try to generate it...')
            if os.path.isfile('depuis_le_abacus.json'):
                print('|-Runtime information: depuis_le_abacus.json found, will use it to generate keywordConversion dictionary')
                keywordConversionJsonFile = jsonRotate(rotateSequence = 'qac', sourceDictionaryFileName='depuis_le_abacus.json')
                with open(keywordConversionJsonFile, 'r') as f:
                    keywordConversion = json.load(f)
            else:
                raise Exception('*ERROR: keywordConversionJsonFile does not exist and depuis_le_abacus.json does not exist')
        else:
            with open(keywordConversionJsonFile, 'r') as f:
                keywordConversion = json.load(f)
    else:
        qeInputScript = {}
        keywordConversion = {}
        raise Exception('*ERROR: mode is not valid')
    discrepantModuleKeywordsExclude = [
        'smearing', 'dftd3_version'
    ]
# INPUT file write from here ======================================================================
    with open(ABACUSInputFile, 'w') as f:
        f.writelines('INPUT_PARAMETERS\n')
        if additionalKeywords["numerical_orbitals"]["use_nao"]:
            f.writelines('basis_type lcao\n')
            print('|-Numerical orbitals: change solver to genelpa')
            qeInputScript["electrons"]["diagonalization"] = "genelpa"
        else:
            f.writelines('basis_type pw\n')
        if len(qeInputScript.keys()) == 0:
            raise Exception('*ERROR: QE input script is empty')

        for section in qeInputScript:

            if section == '__version__':
                if qeInputScript[section] == '>=7.1':
                    qeVersion_7_1 = True
                else:
                    qeVersion_7_1 = False
            elif isNonStandardSectionTitle(section):
                if section != 'CELL_PARAMETERS' and section != 'ATOMIC_POSITIONS' and section != 'K_POINTS' and section != 'ATOMIC_SPECIES':
                    qeNonstandardSections.append(section)
            else:
                for keyword in qeInputScript[section]:
                    if keyword in keywordConversion.keys():
                        if boolUpdateInformation:
                            print('|-Runtime information: original keyword: ' + keyword)
                            print('|-Runtime information: converted keyword: ' + keywordConversion[keyword][0])
                            print('|-Runtime information: pass its value: ' + str(qeInputScript[section][keyword]))
                        convertedKeyword = keywordConversion[keyword][0]
                        value = discrepantKeywordsConversion(keyword, qeInputScript[section][keyword])
                        print('|-Runtime information: comment in depuis_le_abacus.json\n| ' 
                            + ' '*21 + 'Keyword: ' + convertedKeyword + '\n| ' + ' '*21 + 'Comment: ' + keywordConversion[keyword][-1])
                        if convertedKeyword in overwriteKeywords.keys():
                            print('|-Overwrite: will overwrite keyword ' + keyword 
                            + ' with value ' + keywordsWrite(overwriteKeywords[convertedKeyword]))
                            value = overwriteKeywords[convertedKeyword]
                        f.writelines(convertedKeyword + ' ' + keywordsWrite(value) + '\n')
                        # subscript 0 is the position of abacus keyword in the list of values
                    else:
    # Discrapant module conversion ========================================================
                        if keyword == 'occupations':
                            keywordToWrite, valueToWrite = discrepantModuleConversion(qeInputScript, module = 'smearing')
                            for i in range(len(keywordToWrite)):
                                f.writelines(keywordToWrite[i] + ' ' + keywordsWrite(valueToWrite[i]) + '\n')
                        if keyword == 'vdw_corr':
                            keywordToWrite, valueToWrite = discrepantModuleConversion(qeInputScript, module = 'vdw')
                            for i in range(len(keywordToWrite)):
                                f.writelines(keywordToWrite[i] + ' ' + keywordsWrite(valueToWrite[i]) + '\n')
                        print('|-Conversion warning: Keyword ' + keyword + ' not found in conversion table')

        # proceeding non-standard section that unrevelant with STRU and KPT/KLINE files
        for section in qeNonstandardSections:
            if section == 'HUBBARD':
                # well I should say pw is not implemented with dft+u in abacus now...
                # maybe pw is not quite encouraged to be used in abacus? i don't know
                # but i set it here just in case
                if qeVersion_7_1:
                    # because for qe >= 7.1, enabling dft+u calculation is from defining HUBBARD section rather than lda_plus_u
                    # keyword, such keyword is deprecated in version >= 7.1
                    f.writelines('dft_plus_u 1\n')
                else:
                    pass
                
                manifoldCorrection = []
                hubbardUValues = []
                for element in qeInputScript['ATOMIC_SPECIES']['elements']:
                    # abacus also has not implement manifold-by-manifold scheme of dft+u, 
                    # so presently i just use the last one
                    try:
                        manifold = list(qeInputScript['HUBBARD']['on-site'][element].keys())[-1]
                        manifoldCorrection.append(hubbardManifoldToNumber(manifold))
                        hubbardUValues.append(
                            keywordsWrite(
                                qeInputScript['HUBBARD']['on-site'][element][manifold]["U"]
                            )
                            )
                    except KeyError:
                        manifoldCorrection.append('-1')
                        hubbardUValues.append('0.0')
                orbital_corr = ' '.join(manifoldCorrection)
                f.writelines('orbital_corr ' + orbital_corr + '\n')
                hubbard_u = ' '.join(hubbardUValues)
                f.writelines('hubbard_u ' + hubbard_u + '\n')

            elif section == 'SOLVANTS':
                print('|-Conversion warning: a direct conversion from QE-RISM to ABACUS-Implicit solvation model is risky. You should really know what you are doing.')
                print('|-Runtime information: .')
            elif section == 'OCCUPATIONS':
                pass
        for keyword in additionalKeywords:
            if keyword != 'numerical_orbitals':
                f.writelines(keyword + ' ' + keywordsWrite(additionalKeywords[keyword]) + '\n')
            else:
                if additionalKeywords[keyword]["use_nao"]:
                    f.writelines('#numerical_orbitals are added\n')
# STRU file write from here =======================================================================
    with open(ABACUSStructureFile, 'w') as f:
        f.writelines('ATOMIC_SPECIES\n')
        for ityp in range(qeInputScript['system']['ntyp']):
            f.writelines(
                qeInputScript['ATOMIC_SPECIES']['elements'][ityp] 
                + ' ' 
                + keywordsWrite(qeInputScript['ATOMIC_SPECIES']['masses'][ityp])
                + ' ' 
                + keywordsWrite(qeInputScript['ATOMIC_SPECIES']['pseudopotentials'][ityp])
                + '\n')
        # add numerical orbitals
        if additionalKeywords["numerical_orbitals"]["use_nao"]:
            f.writelines('\nNUMERICAL_ORBITAL\n')
            if len(additionalKeywords["numerical_orbitals"]["atom_species"]) == 0:
                print("|-Numerical orbitals: no atom species specified in fp2a.inp, will use atom species from ATOMIC_SPECIES")
                atomSpecies = qeInputScript['ATOMIC_SPECIES']['elements']
            else:
                atomSpecies = additionalKeywords["numerical_orbitals"]["atom_species"]

            naoSelectMode = additionalKeywords["numerical_orbitals"]["select_nao"]
            print("|-Numerical orbitals: select_nao mode: " + naoSelectMode)
            for idxSpecies in range(len(atomSpecies)):
                if naoSelectMode == 'radius_min' or naoSelectMode == 'radius_max':
                    threshold = additionalKeywords["numerical_orbitals"]["cutoff_list_radius"][idxSpecies]
                elif naoSelectMode == 'energy_min' or naoSelectMode == 'energy_max':
                    threshold = additionalKeywords["numerical_orbitals"]["cutoff_list_energy"][idxSpecies]
                elif naoSelectMode == 'energy_min_radius_min' or naoSelectMode == 'energy_min_radius_max' or naoSelectMode == 'energy_max_radius_min' or naoSelectMode == 'energy_max_radius_max':
                    threshold = [
                        additionalKeywords["numerical_orbitals"]["cutoff_list_energy"][idxSpecies], 
                        additionalKeywords["numerical_orbitals"]["cutoff_list_radius"][idxSpecies]
                        ]
                else:
                    raise ValueError("Numerical orbitals: select_nao mode not recognized")
                naoFileName = findNAOByCutoff(
                    element = atomSpecies[idxSpecies],
                    threshold = threshold,
                    zetaNotation = additionalKeywords["numerical_orbitals"]["basis_type"],
                    mode = naoSelectMode
                )
                f.writelines(atomSpecies[idxSpecies] + ' ' + keywordsWrite(naoFileName) + '\n')
        # end of numerical orbitals
        f.writelines('\nLATTICE_CONSTANT\n1.0\n')
        f.writelines('\nLATTICE_VECTORS\n')
        for cell_vector in qeInputScript['CELL_PARAMETERS']['cell']:
            f.writelines(
                keywordsWrite(cell_vector[0]) + ' ' 
                + keywordsWrite(cell_vector[1]) + ' ' 
                + keywordsWrite(cell_vector[2]) + '\n')
        f.writelines('\nATOMIC_POSITIONS\nCartesian\n\n')
        for ityp in range(qeInputScript['system']['ntyp']):
            f.writelines(qeInputScript['ATOMIC_SPECIES']['elements'][ityp] + '\n')
            try:
                if qeInputScript['system']['nspin'] == 1:
                    f.writelines('0.0\n')
                else:
                    try:
                        f.writelines(keywordsWrite(
                            qeInputScript['system']['starting_magnetization({})'.format(str(ityp))]
                            ) + '\n')
                    except KeyError:
                        f.writelines('0.0\n')
            except KeyError:
                f.writelines('0.0\n')
            lineIndices = fromElementFindLinesInXYZ(qeInputScript['ATOMIC_SPECIES']['elements'][ityp], qeInputScript['ATOMIC_POSITIONS']['elements'])
            f.writelines(str(len(lineIndices)) + '\n')
            for lineIndex in lineIndices:
                f.writelines(
                    keywordsWrite(qeInputScript['ATOMIC_POSITIONS']['coordinates'][lineIndex][0])
                    + ' ' + keywordsWrite(qeInputScript['ATOMIC_POSITIONS']['coordinates'][lineIndex][1])
                    + ' ' + keywordsWrite(qeInputScript['ATOMIC_POSITIONS']['coordinates'][lineIndex][2])
                    + ' ')
                f.writelines(
                    keywordsWrite(abs(qeInputScript['ATOMIC_POSITIONS']['constraints'][lineIndex][0]-1))
                    + ' ' + keywordsWrite(abs(qeInputScript['ATOMIC_POSITIONS']['constraints'][lineIndex][1]-1))
                    + ' ' + keywordsWrite(abs(qeInputScript['ATOMIC_POSITIONS']['constraints'][lineIndex][2]-1))
                    + '\n')
            f.writelines('\n')
# KPT  file write from here========================================================================
    with open(ABACUSKpointsFile, 'w', encoding='utf-8') as f:

        f.writelines('K_POINTS\n')
        f.writelines('0\n')
        f.writelines('Gamma\n')
        if qeInputScript['K_POINTS']['mode'] == 'mk':
            f.writelines(
                keywordsWrite(qeInputScript['K_POINTS']['grid'][0]) + ' '
                + keywordsWrite(qeInputScript['K_POINTS']['grid'][1]) + ' '
                + keywordsWrite(qeInputScript['K_POINTS']['grid'][2]) + ' '
                + keywordsWrite(qeInputScript['K_POINTS']['shift'][0]) + ' '
                + keywordsWrite(qeInputScript['K_POINTS']['shift'][1]) + ' '
                + keywordsWrite(qeInputScript['K_POINTS']['shift'][2]) + '\n'
            )
        elif qeInputScript['K_POINTS']['mode'] == 'kpath':
            print('|-K_POINTS warning: In this case you should make sure you already have a KPT file to get converged charge density')
            with open('KLINES', 'w') as f2:
                f2.writelines('K_POINTS\n')
                f2.writelines(str(len(qeInputScript['K_POINTS']['kpts'])) + '\n')
                f2.writelines('Line\n')
                for i, kpt in enumerate(qeInputScript['K_POINTS']['kpts']):
                    f2.writelines(
                        keywordsWrite(kpt[0]) + ' ' 
                        + keywordsWrite(kpt[1]) + ' ' 
                        + keywordsWrite(kpt[2]) + ' ' 
                        + keywordsWrite(qeInputScript['K_POINTS']['w_kpts'][i] )
                        + '\n')
                f2.writelines('\n\n')
    print('*- fp2a (Quantum ESPRESSO-specific version) finished.')

if __name__ == '__main__':
    
    import fp2a_io
    _, overwriteKeywords, additionalKeywords = fp2a_io.readInputScript()

    generateABACUSFiles(
        mode = 'file',
        qeInputScriptJsonFile = "test_qe6_7_in_0.json",
        keywordConversionJsonFile = "depuis_le_abacus_rotated-qac.json",
        qeVersion_7_1 = False,
        overwriteKeywords = overwriteKeywords,
        additionalKeywords = additionalKeywords
    )