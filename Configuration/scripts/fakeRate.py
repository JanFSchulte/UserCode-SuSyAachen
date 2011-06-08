#!/usr/bin/env VERSIONER_PYTHON_PREFER_32_BIT=yes python
'''
Created on 22.12.2010

@author: daniel sprenger, matthias edelhoff
'''

from optparse import OptionParser
from helpers import BetterConfigParser as SafeConfigParser
from copy import copy

from messageLogger import messageLogger as log

import ROOT

def generateFakeRate(cfgFiles):
    from time import asctime
    # ROOT
    ROOT.gROOT.Reset()
    rootContainer = []
    
    log.logInfo("Loading config files %s" % cfgFiles)
    parser = SafeConfigParser()
    parser.read(cfgFiles)
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
                variableToCut = variable
                if parser.has_option("variable:%s"%variable, "abs") and eval(parser.get("variable:%s"%variable, "abs")):
                  variableToCut = "abs(%s)"%variableToCut
                selectionV = "%f <= %s && %s < %f" % (vMin, variableToCut, variableToCut, vMax)
                selectionBin = appendSelection(bin['selection'], selectionV)

                log.logDebug("bin selection: %s" % selectionBin)
                newBin = copy(bin)
                newBin.update({
                               'selection': selectionBin,
                               '%sMin' % variable: vMin,
                               '%sMax' % variable: vMax,
                               })
                newBins.append(newBin)
                if parser.has_option("variable:%s"%variable, "abs") and eval(parser.get("variable:%s"%variable, "abs")):
                    newBin = copy(bin)
                    newBin.update({
                                   'selection': selectionBin,
                                   '%sMin' % variable: -vMax,
                                   '%sMax' % variable: -vMin,
                                   })
                    newBins.append(newBin)  
        bins = newBins

    log.logDebug("Writing to file: %s" % outputFileName)
    file = open(outputFileName, 'w')
    file.write("""#created at %s
#created from  '%s'  

if not "cms" in globals(): import FWCore.ParameterSet.Config as cms


"""%(asctime(), fileName.split()))
    file.close()
    writePSet(outputFileName, parser, tree, bins, "Center")
    writePSet(outputFileName, parser, tree, bins, "Upper")
    writePSet(outputFileName, parser, tree, bins, "Lower")
#    writePSet(outputFileName, parser, tree, bins, "UpperSys")
#    writePSet(outputFileName, parser, tree, bins, "LowerSys")

    return outputFileName

def writePSet(path, parser, trees, bins, name):
    from ROOT import TFile
    import os

    cache = {}
    cacheFilePath = ".fakeRateCache" 
    if os.path.exists(cacheFilePath):
        cacheFile = open(cacheFilePath,"r")
        cache.update( eval(cacheFile.read()) )
        cacheFile.close()
    
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
        loose = []
        tight = []
        weights = []
        for tree in trees:
            inputFileName = tree["path"]
            if not inputFileName in cache:
                cache[inputFileName] = {}
            if not "eventsRunOn" in cache[inputFileName]:
                rootFile = TFile(inputFileName,"read")
                cache[inputFileName]["eventsRunOn"] = rootFile.FindObjectAny("analysis paths").GetBinContent(1)
                rootFile.Close()
            try:
                weights.append( float(eval(parser.get("weights", inputFileName, default = "1.0"),{"eventsRunOn":cache[inputFileName]["eventsRunOn"]}) ))
            except:
                print log.logError("could not convert weight '%s' for '%s'"%(parser.get("weights", inputFileName, default = "1.0"),inputFileName ))
                weights.append(1.0)
            if not selectionLooseBin in cache[inputFileName]:
                cache[inputFileName][selectionLooseBin] = tree["tree"].GetEntries(selectionLooseBin)
            loose.append( int(cache[inputFileName][selectionLooseBin]) )
            
            if not selectionTightBin in cache[inputFileName]:
                cache[inputFileName][selectionTightBin] = tree["tree"].GetEntries(selectionTightBin)
            tight.append( int(cache[inputFileName][selectionTightBin]))
            
        if sum(tight) < 100: log.logWarning("small statistics n= %s for: '%s'"%(tight, selectionTightBin))
        
        (fakeRate, up, low) = getFakeRateWithErrors(tight, loose, weights)
        log.logDebug("%s: %s/%s = %f +%f-%f (weights:%s)" % (bin['selection'], tight, loose, fakeRate,up,low,weights))
        if "Upper" in name:
            fakeRate = up
            
        if "Lower" in name:
            fakeRate=low
            #
          


        # write pset for bin
        file.write("        cms.PSet(\n")
        file.write("            weight = cms.double(%f),\n" % fakeRate)

        for key in bin.keys():
            if (key != "selection"):
                file.write("            %s = cms.double(%f),\n" % (key, bin[key]))
        file.write("        ),\n")
    file.write(")\n")
    file.close()
    
    #write cache
    cacheFile = open(cacheFilePath,"w")
    cacheFile.write(str(cache)+"\n")
    cacheFile.close()

def getFakeRateWithErrors(tight, loose, weights):
    from ROOT import TEfficiency, Double
    from array import array
    from math import sqrt
    assert len(tight) == len(loose) and len(tight) == len(weights), "tight, loose and weights must have same length mave: (%s, %s, %s)"%(len(tight), len(loose), len(weights))
#    print "****"
#    #fakeRate = 1.0 / sum([w*n for w,n in zip(loose, weights)]) * sum([w*n for w,n in zip(tight, weights)])
#    
#    up = Double(-1)
#    low = Double(-1)
#    #uniform prior 
#    fakerate = TEfficiency.Combine(up, low, len(weights), array("i",tight), array("i",loose),1.0,1.0,0.683, array("d",weights),"mode")
#    print "%s\t+ %s\t- %s"%(fakerate, float(low), float(up))
#    
#    wLow = [1./TEfficiency.ClopperPearson(int(loose[i]), int(tight[i]), 0.683, False)**2 for i in range(len(tight))]
#    wUp = [1./TEfficiency.ClopperPearson(int(loose[i]), int(tight[i]), 0.683, True)**2 for i in range(len(tight))]
#    w = [max(low, high) for (low, high) in zip(wLow, wUp)]
#    fakerate = sum([w[i]*tight[i]*1./loose[i] for i in range(len(tight))])/sum([w[i] for i in range(len(tight))])
#    low = sqrt(sum(wLow))
#    up = sqrt(sum(wUp))
#    print "%s\t+ %s\t- %s"%(fakerate, low, up)
    #up = TEfficiency.ClopperPearson(int(loose), int(tight), 0.683, True)
    #low = TEfficiency.ClopperPearson(int(loose), int(tight), 0.683, False)
    fakerate = sum(tight) * 1./sum(loose)
    low = TEfficiency.ClopperPearson(int(sum(loose)), int(sum(tight)), 0.683, False)
    up = TEfficiency.ClopperPearson(int(sum(loose)), int(sum(tight)), 0.683, True)
    #print "%s\t+ %s\t- %s"%(fakerate, low, up)
    
    return (fakerate, up, low)

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
    result = []
    for fileName in fileNames.split():
		tree = ROOT.TChain(treePath)
		tree.Add(fileName)
		if (tree == None):
			log.logError("Could not get tree '%s'\n  from file %s" % (treePath, fileName))
			return None
		else:
			tree.SetDirectory(0)

		result.append({"tree": tree, "path":fileName})
    return result


# entry point
if (__name__ == "__main__"):
    #create option parser
    parser = OptionParser()
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False,
                                  help="talk about everything")
    parser.add_option("-c", "--cfgFile", dest="cfg", action="append", default=[],
                                  help="cfg file(s) to be read")

    (opts, args) = parser.parse_args()
    if (opts.verbose):
        log.outputLevel = 5
    else:
        log.outputLevel = 4

    generateFakeRate(opts.cfg)

