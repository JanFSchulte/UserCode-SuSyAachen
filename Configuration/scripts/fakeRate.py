#!/usr/bin/env VERSIONER_PYTHON_PREFER_32_BIT=yes python

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
    # ROOT
    ROOT.gROOT.Reset()
    rootContainer = []

    log.logInfo("Loading config file %s" % cfgFile)
    parser = SafeConfigParser()
    parser.read(cfgFile)
    variables = parser.sections()
    variables.remove("default")
    log.logDebug("variables: %s" % variables)

    fileName = parser.get('default', 'file')
    outputFileName = parser.get('default', 'outputFile')
    treePath = parser.get('default', 'tree')
    tree = getTreeFromFile(fileName, treePath)
    rootContainer.extend([tree])

    selectionLoose = "ht > 350 && met < 20 && nLept == 1"
    #selectionTight = "ht > 350 && met < 20 && nLept == 1 && tanc > 0.5"
    selectionTight = "ht > 350 && met < 20 && nLept == 1 && pfIso < 0.2"

    log.logInfo("Generating bins")
    bins = [{'selection': ""}]
    for variable in variables:
        binBorders = eval(parser.get(variable, "binBorders"))

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
    file.write("fakes =  cms.VPSet(\n")

    log.logInfo("Looping bins")
    for bin in bins:
        selectionLooseBin = appendSelection(bin['selection'], selectionLoose)
        selectionTightBin = appendSelection(bin['selection'], selectionTight)
        log.logDebug("loose selection: %s" % selectionLooseBin)
        log.logDebug("tight selection: %s" % selectionTightBin)
        loose = tree.GetEntries(selectionLooseBin)
        tight = tree.GetEntries(selectionTightBin)
        fakeRate = 1.0 / loose * tight
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
    return


def appendSelection(selection1, selection2):
    retValue = ""
    if (selection1 == ""):
        retValue = selection2
    elif (selection2 == ""):
        retValue = selection1
    else:
        retValue = "%s && %s" % (selection1, selection2)

    return retValue


def getTreeFromFile(fileName, treePath):
    log.logDebug("Getting tree '%s'\n  from file %s" % (treePath, fileName))

    tree = ROOT.TChain(treePath)
    tree.Add(fileName)

    if (tree == None):
        log.logError("Could not get tree '%s'\n  from file %s" % (treePath, fileName))
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

