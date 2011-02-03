#!/usr/bin/env python
'''
Created on 22.12.2010

@author: daniel sprenger
'''

from optparse import OptionParser
from ConfigParser import SafeConfigParser
from copy import copy

from messageLogger import messageLogger as log

import ROOT


def generateFakeRate(cfgFile):
    from time import asctime
    # ROOT
    ROOT.gROOT.Reset()
    rootContainer = []

    log.logInfo("Loading config file %s" % cfgFile)
    parser = SafeConfigParser()
    parser.read(cfgFile)
    variables = [ i[len("variable:"):] for i in parser.sections() if i.startswith("variable:")]
    log.logDebug("variables: %s" % variables)

    fileName = parser.get('default', 'file')
    outputFileName = parser.get('default', 'outputFile')
    treePath = parser.get('default', 'tree')
    tree = getTreeFromFile(fileName, treePath)
    rootContainer.extend([tree])
    
    log.logInfo("Generating bins")
    bins = [{'selection': ""}]
    for variable in variables:
        binBorders = eval(parser.get("variable:%s"%variable, "binBorders"))

        newBins = []
        for bin in bins:
            for iV in range(0, len(binBorders) - 1):
                vMin = binBorders[iV]
                vMax = binBorders[iV + 1]
                selectionV = "%f <= %s && %s < %f" % (vMin, variable, variable, vMax)
                selectionBin = appendSelection(bin['selection'], selectionV)
                log.logDebug("bin selection: %s" % selectionBin)
                newBin = copy(bin)
                newBin.update({
                               'selection': selectionBin,
                               '%sMin' % variable: vMin,
                               '%sMax' % variable: vMax,
                               })
                newBins.append(newBin)
        bins = newBins

    log.logDebug("Writing to file: %s" % outputFileName)
    file = open(outputFileName, 'w')
    file.write("""#created at %s
#created from  '%s'  

"""%(asctime(), fileName.split()))
    file.close()
    writePSet(outputFileName, parser, tree, bins, "Center")
    writePSet(outputFileName, parser, tree, bins, "Upper")
    writePSet(outputFileName, parser, tree, bins, "Lower")
    return

def writePSet(path, parser, tree, bins, name):
    from ROOT import TEfficiency
    file = open(path, 'a')
    file.write("%s =  cms.VPSet(\n"%(parser.get("default","outputPSet", raw=True)%name))

    selectionLoose = parser.get("selections","selectionLoose") 
    selectionTight = parser.get("selections","selectionTight")
   
    log.logInfo("Looping bins")
    for bin in bins:
        selectionLooseBin = appendSelection(bin['selection'], selectionLoose)
        selectionTightBin = appendSelection(bin['selection'], selectionTight)
        log.logDebug("loose selection: %s" % selectionLooseBin)
        log.logDebug("tight selection: %s" % selectionTightBin)
        loose = tree.GetEntries(selectionLooseBin)
        tight = tree.GetEntries(selectionTightBin)
        if tight < 100: log.logWarning("small statistics n= %s for: '%s'"%(tight, bin["selection"]))
        fakeRate = 1.0 / loose * tight
        if(name == "Upper"):
          fakeRate = TEfficiency.ClopperPearson(loose, tight, 0.683, True)
        if(name == "Lower"):
          fakeRate = TEfficiency.ClopperPearson(loose, tight, 0.683, False)
          
        log.logDebug("%s: %d/%d = %f" % (bin['selection'], tight, loose, fakeRate))

        # write pset for bin
        file.write("        cms.PSet(\n")
        file.write("            weight = cms.double(%f),\n" % fakeRate)

        for key in bin.keys():
            if (key != "selection"):
                file.write("            %s = cms.double(%f),\n" % (key, bin[key]))
        file.write("        ),\n")
    file.write(")\n")
    file.close()


def appendSelection(selection1, selection2):
    retValue = ""
    if (selection1 == ""):
        retValue = selection2
    elif (selection2 == ""):
        retValue = selection1
    else:
        retValue = "%s && %s" % (selection1, selection2)

    return retValue


def getTreeFromFile(fileNames, treePath):
    log.logDebug("Getting tree '%s'\n  from file %s" % (treePath, fileNames.split()))

    tree = ROOT.TChain(treePath)
    for fileName in fileNames.split():
      tree.Add(fileName)

    if (tree == None):
        log.logError("Could not get tree '%s'\n  from file %s" % (treePath, fileNames.split()))
        return None
    else:
        tree.SetDirectory(0)
        return tree


# entry point
if (__name__ == "__main__"):
    #create option parser
    parser = OptionParser()
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False,
                                  help="talk about everything")
    parser.add_option("-c", "--cfgFile", dest="cfg", action="store", type="string", default="fakeRate.ini",
                                  help="cfg file to be read")

    (opts, args) = parser.parse_args()
    if (opts.verbose):
        log.outputLevel = 5
    else:
        log.outputLevel = 4

    generateFakeRate(opts.cfg)

