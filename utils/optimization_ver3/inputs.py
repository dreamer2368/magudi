from warnings import warn
import yaml

verbose = True

class InputParser:
    dict_ = None
    filename = ""

    def __init__(self, filename):
        with open(filename) as c:
            self.dict_ = yaml.safe_load(c)
        self.filename = filename
        return

    def getInput(self, keys, fallback=None, datatype=None):
        keyString = ""
        for key_ in keys:
            keyString += key_ + "/"

        val = self.dict_
        for key in keys:
            if key in val:
                val = val[key]
            elif (fallback != None):
                return fallback
            else:
                raise RuntimeError("%s does not exist in the input file %s!" % (keyString, self.filename))

        if (fallback != None):
            if (type(val) != type(fallback)):
                raise RuntimeError("%s does not match the type with the fallback value %s!" % (str(type(val)), str(type(fallback))))
        elif (datatype != None):
            if (type(val) != datatype):
                raise RuntimeError("%s does not match the specified datatype %s!" % (str(type(val)), str(datatype)))
        else:
            from inputs import verbose
            if verbose: warn("InputParser Warning: datatype is not checked.\n key: %s\n value type: %s" % (keys, type(val)))
        return val

# get a dict with {key: val} from a list of dicts
# NOTE: it returns only the first item in the list,
#       even if the list has more than one dict with {key: val}.
def getDictFromList(list_, inputDict):
    dict_ = None
    for item in list_:
        isDict = True
        for key, val in inputDict.items():
            if key not in item:
                isDict = False
                break
            if (item[key] != val):
                isDict = False
                break
        if (isDict):
            dict_ = item
            break
    if (dict_ == None):
        raise RuntimeError('Given list does not have a dict with {%s: %s}!' % (key, val))
    return dict_
