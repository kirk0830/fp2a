import re
def keywordsWrite(keyword):

    # datatype converter
    if type(keyword) == str:
        return keyword
    elif type(keyword) == bool:
        return str(int(keyword))
    else:
        try:
            return str(keyword)
        except ValueError:
            raise ValueError('ERROR: ' + keyword + ' is not a valid value for keyword ' + keyword)

def keywordsRead(keyword):
    """
    keywordPostProcessing is for post-processing the keyword value based on the type of the keyword value.
    For a number which should be float but the format of read-in is string, this function will convert it to float.\n
    @param keyword: the keyword value read-in from the input script\n
    @return: the post-processed keyword value\n
    """
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
                scientificNotation = re.split('([eEdD])', keyword)
                try:
                    return float(scientificNotation[0]) * 10 ** int(scientificNotation[2])
                except ValueError:
                    return keyword
        else:
            try:
                return int(keyword)
            except ValueError:
                return keyword
                