import re
def readInputScript(fileName = 'fp2a.inp'):
    # read input script
    f = open(fileName, 'r')
    lines = f.readlines()
    f.close()
    # parse input script
    basicConversion = {}
    additionalKeywords = {}
    overwriteKeywords = {}
    presentSection = 'nothing'
    for line in lines:

        cleanLine = re.sub('\s|\t|\n', '', line.strip())
        sharpIndex = cleanLine.find('#')
        if sharpIndex != -1:
            cleanLine = cleanLine[:sharpIndex]
        cleanLine = cleanLine.replace(' ', '')
        #print('line starts->|'+cleanLine+'|<-line ends')
        if len(cleanLine) == 0:
            continue
        else:
            words = re.split(',|=', cleanLine)
            if words[0] == '&basic_conversion':
                presentSection = 'basicConversion'
            elif words[0] == '&additional_keywords':
                presentSection = 'additionalKeywords'
            elif words[0] == '&overwrite_keywords':
                presentSection = 'overwriteKeywords'
            elif words[0] == '/':
                continue
            else:
                if presentSection == 'basicConversion':
                    if words[0] != 'script_name':
                        basicConversion[words[0]] = words[1]
                    else:
                        basicConversion[words[0]] = words[1:]
                elif presentSection == 'additionalKeywords':
                    additionalKeywords[words[0]] = words[1]
                elif presentSection == 'overwriteKeywords':
                    overwriteKeywords[words[0]] = words[1]

    return basicConversion, overwriteKeywords, additionalKeywords

if __name__ == '__main__':
    basicConversion, overwriteKeywords, additionalKeywords = readInputScript()
    print(basicConversion)
    print(overwriteKeywords)
    print(additionalKeywords)