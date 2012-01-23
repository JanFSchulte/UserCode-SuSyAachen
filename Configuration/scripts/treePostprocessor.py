#!/usr/bin/env python
'''
Created on 26.05.2011

@author: heron
'''

class TreeProcessor:
    def __init__(self, config, name, isMC = False):
        self.config = config
        self.section = "treeProcessor:%s"%name
        self.isMC = isMC
        self.nThEntry = 0
        self.nEntries = None
        
    def prepareSrc(self, tree, object, processors):
        self.nEntries = tree.GetEntries()
        #print "%s: %d" % (self.section, self.nEntries)
        
    
    def prepareDest(self, tree, object):
        pass

    def processEvent(self, event, object):
        import sys
        #self.nThEntry +=1
        #if self.nThEntry%(self.nEntries*0.001) < 1:
        #    #print "%d / %d" % (self.nThEntry, self.nEntries)
        #    sys.stdout.write("\r%.1f%%" %(self.nThEntry*1./(self.nEntries*0.01)))    
        #    sys.stdout.flush()  
        return True
    
class SimpleSelector(TreeProcessor):
    def __init__(self, config, name, isMC = False):
        TreeProcessor.__init__(self, config, name, isMC)
        
    def processEvent(self, event, object):
        from math import fabs
        expression = self.getExpression(object)
        expression = expression.replace("&&","and")
        expression = expression.replace("&","and")
        expression = expression.replace("||","or")
        expression = expression.replace("|","or")        
        evalGlobal = {"abs":fabs}
        for i in [i.GetName() for i in event.GetListOfBranches()]:
            evalGlobal[i] = getattr(event,i)
        return eval(expression, evalGlobal)

    def getExpression(self, object):
        #print self.section
        #print "%sExpression"%object
        #print self.config.get(self.section,"%sExpression"%object)
        return self.config.get(self.section,"%sExpression"%object)
        
    
class SimpleWeighter(TreeProcessor):
    def __init__(self, config, name, isMC = False):
        TreeProcessor.__init__(self, config, name, isMC)
        self.weight = {}
        
    def prepareSrc(self, src, object, processors):
        #src.SetBranchStatus("weight", 0)
        pass
        
    def prepareDest(self, dest, object):
        from array import array
        self.weight[object] = array("f",[1.0])
        dest.Branch("weight2",self.weight[object],"weight2/F")

    def processEvent(self, event, object):
        self.weight[object][0] = event.weight*(event.pt1+event.pt2)
        return True

class NoFakeWeighter(TreeProcessor):
    def __init__(self, config, name, isMC = False):
        TreeProcessor.__init__(self, config, name, isMC)
        self.branchName = config.get(self.section, "branchName")
        self.noFakeValue = eval(config.get(self.section, "weight"))
        self.weight = {} 

    def prepareSrc(self, src, object, allProcessors):
        TreeProcessor.prepareSrc(self, src, object, allProcessors)
        src.SetBranchStatus(self.branchName, 0)
        
    def prepareDest(self, dest, object):
        from array import array
        TreeProcessor.prepareDest(self, dest, object)
        self.weight[object] = array("f",[self.noFakeValue])
        dest.Branch(self.branchName,self.weight[object],"%s/F"%self.branchName)
        
    #def processEvent(self, event, object):
    #    TreeProcessor.processEvent(self, event, object)
    #    return True


class FakeWeighter(TreeProcessor):
    def __init__(self, config, name, isMC = False):
        from helpers import parsePSet
        TreeProcessor.__init__(self, config, name, isMC)
        frPath = {}
        self.cuts = {}
        self.pSet = {}
        for optionName in config.options(self.section):
            if optionName.startswith("cut_"):
                name = optionName.split("cut_")[1]
                frPath[name] = config.get(self.section, "fakePSet_%s"%name)
                self.cuts[name] =  config.get(self.section, "cut_%s"%name)
                assert not frPath[name] == None, "could not find path 'fakePSet_%s' for cut 'cut_%s' in section '%s'"%(name, name, self.section) 
        if config.has_option(self.section, "fakePSet"):
            frPath["data"] = config.get(self.section, "fakePSet")
            self.cuts["data"] = "True"    
            frPath["MC"] = frPath["data"]
        self.branchName = config.get(self.section, "branchName")
        for name, path in frPath.items():            
            self.pSet[name] = parsePSet(path)
        #self.pSetMC = parsePSet(mcFrPath)
        self.ptMax = 74.9
        self.weight = {}
        
    def prepareSrc(self, src, object, allProcessors):
        TreeProcessor.prepareSrc(self, src, object, allProcessors)
        src.SetBranchStatus(self.branchName, 0)
        
    def prepareDest(self, dest, object):
        from array import array
        TreeProcessor.prepareDest(self, dest, object)
        self.weight[object] = array("f",[-1.0])
        dest.Branch(self.branchName, self.weight[object],"%s/F"%self.branchName)
        
    def processEvent(self, event, object):
        TreeProcessor.processEvent(self, event, object)
        frBins = {}
        for (binName, varName) in zip(self.config.get(self.section,"binNames").split(),self.config.get(self.section,"varNames").split()):
            frBins[binName] = getattr(event,varName)
        for binName in frBins:
            if "pt" in binName and frBins[binName] > self.ptMax:
                frBins[binName] = self.ptMax
        idExpr = self.config.get(self.section,"idExpr")
        evalGlobals = {"id1":event.id1,
                       "id2":event.id2,
                       "runNr": event.runNr,
                       "isMC": self.isMC}
        self.weight[object][0] = -1

        usedCut = None
        for cutName, cut in self.cuts.items():
            if eval(cut, evalGlobals):
                assert usedCut == None, "event already processed by '%s: %s', also passed '%s: %s'"%(usedCut, self.cuts[usedCut],
                                                                                                              cutName,self.cuts[cutName])
                usedCut = cutName
        if usedCut == None:
            usedCut = "data"
            
        if eval(idExpr, evalGlobals):
            frType = "center"
            frType = self.config.get(self.section,"frType")
            pSet = self.pSet[usedCut]
            assert frType in pSet, "in '%s' could not find fakerate '%s' in '%s'"%(self.section, frType, pSet)
                
            f = pSet[frType].fakeRate(frBins)
            self.weight[object][0] = f/(1-f)
        return True

class Weighter(TreeProcessor):
    def __init__(self, config, name, isMC = False):
        from helpers import parsePSet
        TreeProcessor.__init__(self, config, name, isMC)
        self.weightPSets = {}
        for pSetName in ["weight1PSet","weight2PSet","weight1PSetMC","weight2PSetMC"]:            
            self.weightPSets[pSetName] = None
            if config.has_option(self.section, pSetName):
                pSetPath = config.get(self.section, pSetName)
                self.weightPSets[pSetName] = parsePSet( pSetPath )
        self.branchName = config.get(self.section, "branchName")
        self.formula = config.get(self.section, "formula")
        self.binMax = {}
        for optName in self.config.options(self.section):
            if optName.endswith("Max"):
                self.binMax[optName[:-3]] = eval(self.config.get(self.section, optName))
        self.weight = {}

    def prepareSrc(self, src, object, allProcessors):
        TreeProcessor.prepareSrc(self, src, object, allProcessors)
        src.SetBranchStatus(self.branchName, 0)
        
    def prepareDest(self, dest, object):
        from array import array
        TreeProcessor.prepareDest(self, dest, object)
        self.weight[object] = array("f",[-1.0])
        dest.Branch(self.branchName, self.weight[object],"%s/F"%self.branchName)
        
    def processEvent(self, event, object):
        TreeProcessor.processEvent(self, event, object)
        frBins = {}
        for (binName, varName) in zip(self.config.get(self.section,"binNames").split(),self.config.get(self.section,"varNames").split()):
            frBins[binName] = getattr(event,varName)
        print frBins
        for binName in frBins:
            if binName in self.binMax and frBins[binName] > self.binMax[binName]:
                frBins[binName] = self.binMax[binName]
        idExpr = self.config.get(self.section,"idExpr")
        evalGlobals = {"id1":event.id1,
                       "id2":event.id2}
        self.weight[object][0] = -1
        if eval(idExpr,evalGlobals):
            frType = "center"
            frType = self.config.get(self.section,"frType")
            pSet = self.pSet
            if self.isMC: pSet = self.pSetMC
            assert frType in pSet, "in '%s' could not find fakerate '%s' in '%s'"%(self.section, frType, pSet)
        
            f = pSet[frType].fakeRate(frBins)
            self.weight[object][0] = f/(1-f)
        return True

    
class OverlapRemover(TreeProcessor):
    def __init__(self, config, name, isMC = False):
        from os.path import exists as pathExists
        from os import makedirs
        TreeProcessor.__init__(self, config, name, isMC)
        self.keepEvents = {}
        self.rejected = {}
        self.listPath = self.config.get(self.section,"listPath")
        if not pathExists(self.listPath):
            makedirs(self.listPath)
    
    def prepareSrc(self, src, object, allProcessors):
        TreeProcessor.prepareSrc(self, src, object, allProcessors)
        endOfLine = 10
        for ev in src:
            if (endOfLine < 1):
                pass
                #continue
            endOfLine -= 1
            processingResults = {}
            processors = self.config.get(self.section,"%sProcessors"%object).split()
            filter = " and ".join(processors)
            if self.config.has_option(self.section,"%sFilter"%object):
                filter = self.config.get(self.section,"%sFilter"%object)
            for processorName in processors:
                processingResults[processorName] = allProcessors[processorName].processEvent(ev, object)
            if eval(filter, processingResults):
                self._processEvent(ev, object)
        self._writeEventList("%s/selected.eventList"%(self.listPath), self.keepEvents.keys())
        for rejectedObject, rejectedList in self.rejected.iteritems():
            self._writeEventList("%s/rejected.%s.eventList"%(self.listPath, rejectedObject), rejectedList)
                            
    def _processEvent(self, ev, object):
        rejectedObject = None
        fingerPrint = (ev.runNr, ev.lumiSec, ev.eventNr) 
        if not fingerPrint in self.keepEvents:
            self.keepEvents[fingerPrint] = (object, ev.pt1+ev.pt2)
        elif self.keepEvents[fingerPrint][1] < ev.pt1+ev.pt2 or ("Tau" in self.keepEvents[fingerPrint][0] and not "Tau" in object):
            rejectedObject = self.keepEvents[fingerPrint][0]
            self.keepEvents[fingerPrint] = (object, ev.pt1+ev.pt2)                
        else:
            rejectedObject = object
        if rejectedObject:
            if not rejectedObject in self.rejected:
                self.rejected[rejectedObject] = []
            self.rejected[rejectedObject].append(fingerPrint)     
            
    def _writeEventList(self, path, list):
            file = open(path,"w")
            file.write("\n".join([":".join(["%s"%j for j in i]) for i in list]))
            file.close()
    
    def processEvent(self, event, object):
        TreeProcessor.processEvent(self, event, object)
        fingerPrint = (event.runNr, event.lumiSec, event.eventNr)        
        #if fingerPrint in self.keepEvents and not object == self.keepEvents[fingerPrint][0]:
        #    print "skipping", fingerPrint, object, event.pt1+event.pt2,   self.keepEvents[fingerPrint]
        value = fingerPrint in self.keepEvents and object == self.keepEvents[fingerPrint][0]
        value2 = True
        if (value):
            # checks if correct sumPt is matched
            # necessary for same-flavour events
            value2 = event.pt1+event.pt2 == self.keepEvents[fingerPrint][1]
        
        return (value and value2)

class TreeProducer:
    def __init__(self, config, inputPaths, name = None, isMC = False):
        from os.path import split as splitPath
        from os.path import exists as pathExists
        from os import makedirs
        from sys import modules
        if name == None:
            assert len(inputPaths) == 1, "can only determine names for single paths automatically. Got '%s'"%inputPaths
            name = splitPath(inputPaths[0])[1].split(".")[2]
        self.config = config
        self.name = name
        self.isMC = isMC 
        self.tasks = list(set([splitPath(i)[1].split(".")[1] for i in inputPaths]))
        self.flags = list(set([splitPath(i)[1].split(".")[0] for i in inputPaths]))
        self.inputPaths = inputPaths
        self.counterSum = None
        self.outPath = config.get("general","outPath")
        self.outTask = ""
        if self.config.has_option("general","outTask"):
            self.outTask = self.config.get("general","outTask")
        if not pathExists(self.outPath):
            makedirs(self.outPath)
        self.treeProcessors = {}
        for section in self.config.sections():
            if section.startswith("treeProcessor:"):
                processorName = section.split("treeProcessor:")[1]
                processorType = self.config.get(section,"type")
                #warning black magic ahead :P
                self.treeProcessors[processorName] = getattr(modules[globals()["__name__"]],processorType)(self.config, processorName, isMC = self.isMC)
    
    def __repr__(self):
        from pprint import pformat
        result = "TreeProducer: %s %s \n"%(self.name, "(MC)" if self.isMC else "(Data)")
        for section in self.config.sections():            
            trees = self.getTrees(section)
            if trees:
                result += "\n%s:\n"%section 
                result += pformat(trees)        
        return result
        
    def produce(self):
        from ROOT import TFile, TChain, TH1
        from os.path import exists as pathExists
        from os.path import split as splitPath
        outFilePath = "%s/%s.%s.%s.root"%(self.outPath, "".join(self.flags), "processed", self.name)
        if pathExists(outFilePath):
            print "skipping", outFilePath
            return
        outFile = TFile("%s/%s.%s.%s.root"%(self.outPath, "".join(self.flags), "processed", self.name),"recreate")
        
        for section in self.config.sections():
            treeTuple= self.getTrees(section)                           
            if not treeTuple == False:
                (trees, treeName, subDirName) = treeTuple
                outDir = None
                srcTree = {} 
                for object in trees:
                    srcTree[object] = TChain("%s%s"%(object, treeName))
                    processors = self.config.get(section,"%sProcessors"%object).split()
                    filter = " and ".join(processors)
                    if self.config.has_option(section,"%sFilter"%object):
                        filter = self.config.get(section,"%sFilter"%object)
                    for treePath in trees[object]:
                        #srcFile = TFile(filePath,"r")
                        #srcTree = srcFile.Get(treePath)
                        filePath = "%s.root"%treePath.split(".root")[0]
                        inFile = TFile(filePath,"READ")
                        if not self.counterSum:
                            dirName = "%sCounters" % section.split("dileptonTree:")[1]
                            if not self.outTask == "":
                                dirName = "%sCounters"%self.outTask
                            outFile.mkdir(dirName)
                            outFile.cd(dirName)
                            task = None                            
                            for t in self.tasks:
                                if ".%s."%t in splitPath(filePath)[1]:
                                    assert task == None, "unable to disambiguate tasks '%s' matches both '%s' and '%s'"(filePath, task, t)
                                    task = t                        
                            self.counterSum = inFile.Get("%sCounters/analysis paths"%task).Clone()
                        else:
                            pass
                            #need to cope with different lumis :( 
                            #h = inFile.Get("%sCounters/analysis paths"%task)
                            #print inFile, "%sCounters/analysis paths"%task, h
                            #self.counterSum.Add( h,1. )
                        inFile.Close()
                        srcTree[object].Add(treePath)
                        print "adding", treePath
                    srcTree[object].SetBranchStatus("*", 1)
                    for processorName in processors:
                        if (self.treeProcessors[processorName].__class__.__name__ == SimpleSelector.__name__ and not self.config.has_option(section,"%sFilter"%object)):
                            print "Requirements met, applying simple selection boosting =)"
                            expression = self.treeProcessors[processorName].getExpression(object)
                            print "Cutting tree down to: '%s'" % (expression)
                            srcTree[object] = srcTree[object].CopyTree(expression)
                        self.treeProcessors[processorName].prepareSrc(srcTree[object], object, self.treeProcessors)
                for object in trees:
                    print object
                    processors = self.config.get(section,"%sProcessors"%object).split()
                    filter = " and ".join(processors)
                    if self.config.has_option(section,"%sFilter"%object):
                        filter = self.config.get(section,"%sFilter"%object)
                    
                    if not outDir:
                        outDir = outFile.mkdir(subDirName)
                    outFile.cd(subDirName)
                    destTree = srcTree[object].CloneTree(0)
                    #print processors
                    for processorName in processors:
                        self.treeProcessors[processorName].prepareDest(destTree, object)
                        print "%s: %d" % (str(processorName), self.treeProcessors[processorName].nEntries)
                    endOfLine = 1000
                    for i in srcTree[object]:
                        if endOfLine < 1:
                            pass
                            #continue
                        endOfLine -= 1
                        processingResults = {}
                        for processorName in processors:
                            processingResults[processorName] = self.treeProcessors[processorName].processEvent(srcTree[object], object)
                        if eval(filter, processingResults):
                            destTree.Fill()
                    #srcFile.Close()
                    outFile.Write()
                #from pprint import pprint
                #pprint( trees)
        
        outFile.Close()
    
    def openTrees(self):
        from ROOT import TChain
        for section in self.config.sections():
            result = {}
            treeTuple= self.getTrees(section)                           
            if not treeTuple == False:
                (trees, treeName, subDirName) = treeTuple
                for channel, paths in trees.items():
                    result[channel] = TChain(paths[0].split(".root/")[1])
                    for path in paths:
                        result[channel].Add(path)
                yield section, result
        
        
    def getTrees(self, section):
        trees = None
        if section.startswith("dileptonTree:"):
            treeProducerName =self.config.get(section,"treeProducerName")
            trees = self._getDileptonTrees(section)
            treeName = "DileptonTree"
            subDirName = "%s%s%s"%(self.outTask, section.split("dileptonTree:")[1], treeProducerName)
        if section.startswith("isoTree:"):
            treeProducerName =self.config.get(section,"treeProducerName")
            trees = self._getIsoTrees(section)
            treeName = "Iso"
            subDirName = "%s%s"%(section.split("isoTree:")[1],treeProducerName)
        if trees == None:
            return False
        else:
            return trees, treeName, subDirName
    
    def _getIsoTrees(self, section):
        from re import match
        from os.path import split as splitPath
        result = {}
        datasetPaths = []
        treeProducerName =self.config.get(section,"treeProducerName")
        datasetExpression = self.config.get(section, "Dataset")
        datasetSelection = self.config.get(section, "Selection")
        for path in self.inputPaths:
            if not match(datasetExpression, splitPath(path)[1].split(".")[2]) == None:
                datasetPaths.append(path)
        if datasetPaths == []:
            if self.config.has_option(section, "otherSelection"):  
                datasetSelection = self.config.get(section, "otherSelection")
            datasetPaths = self.inputPaths
        result[""] = []
        for path in datasetPaths:
            task = None
            for t in self.tasks:
                if ".%s."%t in splitPath(path)[1]:
                    assert task == None, "unable to disambiguate tasks '%s' matches both '%s' and '%s'"(path, task, t)
                    task = t
            result[""].append( "%s/%s%s%s/Trees/Iso"%(path,task,datasetSelection,treeProducerName))
        return result
        
        
    def _getDileptonTrees(self, section):
        from re import match
        from os.path import split as splitPath
        result = {}
        objects = self.config.get(section,"objects").split()
        treeProducerName =self.config.get(section,"treeProducerName")
        for object in objects:
            result[object] = []
            datasetExpressions = self.config.get(section, "%sDataset"%object).split()
            datasetSelections = self.config.get(section, "%sSelection"%object).split()
            if datasetSelections == []: datasetSelections = [""] *len(datasetExpressions)
            objectProducerNames = [treeProducerName]*len(datasetExpressions)
            if self.config.has_option(section, "%sTreeProducerName"%object):
                objectProducerNames = self.config.get(section, "%sTreeProducerName"%object).split()
                if objectProducerNames == []: objectProducerNames = [""] *len(datasetExpressions)
            assert len(datasetExpressions) == len(datasetSelections), "length missmatch in selections '%s' <> '%s' "%(datasetExpressions, datasetSelections) 
            assert len(datasetExpressions) == len(objectProducerNames), "length missmatch in producer names'%s' <> '%s' "%(datasetExpressions, objectProducerNames)
            for path in self.inputPaths:
                datasetSelection = None
                objectProducerName = None
                for i in range(len(datasetExpressions)):
                    datasetExpression = datasetExpressions[i]                                    
                    if not match(datasetExpression, splitPath(path)[1].split(".")[2]) == None:
                        datasetSelection = datasetSelections[i]
                        objectProducerName = objectProducerNames[i]
                if datasetSelection == None:
                    if self.config.has_option(section, "%sMCSelection"%object) and self.isMC:
                        datasetSelection = self.config.get(section, "%sMCSelection"%object)
                        objectProducerName = treeProducerName
                        print "using MC selection", datasetSelection 
#                    elif self.config.has_option(section, "otherSelection"):
#                        datasetSelection = self.config.get(section, "otherSelection")
#                        objectProducerName = treeProducerName    
                task = None
                for t in self.tasks:
                    if ".%s."%t in splitPath(path)[1]:
                        assert task == None, "unable to disambiguate tasks '%s' matches both '%s' and '%s'"(path, task, t)
                        task = t
                if not datasetSelection == None:
                    result[object].append( "%s/%s%s%s/%sDileptonTree"%(path,task,datasetSelection,objectProducerName,object))
            assert len(result[object]) > 0, "'otherSelection' has been depricated, please use MCSelection instead! (for %s matching '%s' in '%s')"%(object, 
                                                                                                                                                    self.config.get(section, "%sDataset"%object), 
                                                                                                                                                    [splitPath(path)[1].split(".")[2] for path in self.inputPaths])
                                                                                                                                                    
                
        return result
        
def getProducers(config, rawPaths):
    from glob import glob
    from re import match
    from os.path import split as splitPath
    from os.path import isdir
    mcExpression = config.get("general","MCDatasets")
    result = []
    dataPaths = []
    paths = []
    for path in rawPaths:
        if isdir(path): 
            path = "%s/*.root"%path
            paths.extend( glob(path) )
        else:
            paths.append(path)
    for inputPath in paths:
        if match(mcExpression, splitPath(inputPath)[1]) ==None:
            dataPaths.append(inputPath)
        else:   
            result.append(TreeProducer(config, [inputPath], isMC=True))
    if len(dataPaths) > 0:
        result.append(TreeProducer(config, dataPaths, "MergedData", isMC=False))
    return result

    
    
def main(argv = None):
    import sys
    #from ConfigParser import ConfigParser
    from helpers import BetterConfigParser as ConfigParser
    from optparse import OptionParser
    if argv == None:
        argv = sys.argv[1:]
    parser = OptionParser()
    parser.add_option("-C", "--Config", dest="Config", action="append", default=[],
                          help="Main configuration file. Can be given multiple times in case of split configurations. Default is Input/default.ini")
    parser.add_option("-p", "--path", dest="paths", action="append", default=[],
                          help="paths to consider (wildcards allowed, but but them in ''!)")
    (opts, args) = parser.parse_args(argv)

    if opts.Config == []:
        opts.Config = [ "default.ini" ]
    
    config = ConfigParser()
    config.read(opts.Config)

    basePaths = opts.paths
    if opts.paths == []:
        basePaths = [config.get("general","basePath")]
    
    producers = []
     
    producers.extend( getProducers(config, basePaths))
    for p in producers:
        p.produce()
    
if __name__ == '__main__':
    main()

import unittest    

class plotTest(unittest.TestCase):
    configString ="""
[general]
tasks = pfDiLeptonCuts triggerStudies
basePath = /Users/heron/Documents/superSymmetry/results/diLeptonTaus/diLeptons42x/largesusy0453v14/input/
#/Users/heron/Documents/superSymmetry/results/diLeptonTaus/diLeptons41x/susy0447v5/input
MCDatasets = .*_Summer11
outPath = processedTrees

[dileptonTree NoCuts]
treeProducerName = TaNCTrees
objects = EE EMu MuMu ETau MuTau TauTau
EEDataset = SingleElectron.* 
EESelection = 
EEProcessors = htSelector ptSumWeighter overlap
EEFilter = True
EMuDataset = SingleMu_.*
EMuSelection = 
EMuProcessors = htSelector ptSumWeighter overlap
EMuFilter = htSelector 
MuMuDataset = SingleMu_.*
MuMuSelection = 
MuMuProcessors = htSelector overlap
ETauDataset = TauPlusX_.*
ETauSelection = HLTETau
ETauProcessors = htSelector overlap
MuTauDataset = TauPlusX_.*
MuTauSelection = HLTMuTau
MuTauProcessors = htSelector overlap
TauTauDataset = TauPlusX_.*
TauTauSelection = HLTTauTauMHT
TauTauProcessors = htSelector overlap
OtherSelection =

[dileptonTree NoCuts]
treeProducerName = HpsACTrees
objects = EE EMu MuMu ETau MuTau TauTau
EEDataset = SingleElectron.* 
EESelection = 
EEProcessors = overlap
EMuDataset = SingleMu_.*
EMuSelection = 
EMuProcessors = overlap 
MuMuDataset = SingleMu_.*
MuMuSelection = 
MuMuProcessors = overlap
ETauDataset = TauPlusX_.*
ETauSelection = HLTETau
ETauProcessors = overlap
MuTauDataset = TauPlusX_.*
MuTauSelection = HLTMuTau
MuTauProcessors = overlap
TauTauDataset = TauPlusX_.*
TauTauSelection = HLTTauTauHT
TauTauProcessors = overlap
OtherSelection = 
 
[dileptonTree:HighPt]
treeProducerName = HpsACTrees 
objects = EE EMu MuMu ETau MuTau TauTau
EETreeProducerName = FinalTrees 
EEDataset = Data_Run2011.*
EESelection = 
EEProcessors = overlap
#EEFilter = singleLooseNotTightSelector or overlapHighPt
EMuTreeProducerName = FinalTrees 
EMuDataset = Data_Run2011.*
EMuSelection =  
EMuProcessors = overlap 
#EMuFilter = singleLooseNotTightSelector or overlapHighPt
MuMuTreeProducerName = FinalTrees 
MuMuDataset = Data_Run2011.*
MuMuSelection = 
MuMuProcessors = overlap 
#MuMuFilter = singleLooseNotTightSelector or overlapHighPt
ETauDataset = TauPlusX_.*
#ETauSelection = 
ETauSelection = HLTETau
ETauProcessors = overlap 
#ETauFilter = singleLooseNotTightSelector or overlapHighPt
MuTauDataset = TauPlusX_.*
#MuTauSelection = 
MuTauSelection = HLTMuTau
MuTauProcessors = overlap 
#MuTauFilter = singleLooseNotTightSelector or overlapHighPt
TauTauDataset = TauPlusX_.*
TauTauMCSelection = 
#TauTauSelection = 
TauTauSelection = HLTTauTauMHT
TauTauProcessors = overlap 
#TauTauFilter = singleLooseNotTightSelector or overlapHighPt
#otherSelection = 

 

[isoTree NoCuts]
treeProducerName = TnPTaNCTauTrees
Dataset = TauPlusX
Selection = 
Processors = tauFakeWeights
PtherSelection = 

[treeProcessor:htSelector]
type = SimpleSelector

[treeProcessor:lowPtSelector]
type = SimpleSelector
EEExpression = pt1 > 20 && pt2 > 20 && ht > 250 
EMuExpression = pt1 > 20 && pt2 > 20 && ht > 250 
MuMuExpression = pt1 > 20 && pt2 > 20 && ht > 250 
ETauExpression = pt1 > 20 && pt2 > 15 && ht > 250 
MuTauExpression = pt1 > 20 && pt2 > 15 && ht > 250 
TauTauExpression = pt1 > 20 && pt2 > 15 && ht > 250 


[treeProcessor:ptSumWeighter]
type = SimpleWeighter

[treeProcessor:tauFakeWeights]
type = FakeWeighter
fakePSet = /Users/heron/Documents/superSymmetry/results/diLeptonTaus/diLeptons41x/susy0447v5/tauDataHT_cff.py
idExpr = id2 < 0.
varNames = pt2 eta2
binNames = pt eta
branchName = fakeWeight2


[treeProcessor:tauWeights]
type = Weighter
weight1PSet = /Users/heron/Documents/superSymmetry/results/diLeptonTaus/diLeptons41x/susy0447v5/tauDataHT_cff.py
weight2PSet = /Users/heron/Documents/superSymmetry/results/diLeptonTaus/diLeptons41x/susy0447v5/tauDataHT_cff.py
weight1PSetMC = /Users/heron/Documents/superSymmetry/results/diLeptonTaus/diLeptons42x/largeSusy0452v1/hpsTauQCDMC_cff.py
weight2PSetMC = /Users/heron/Documents/superSymmetry/results/diLeptonTaus/diLeptons42x/largeSusy0452v1/hpsTauQCDMC_cff.py

idExpr = id2 < 0.
varNames = pt2 eta2
binNames = pt eta
ptMax = 75.
branchName = fakeWeight2
formula = w1 * w2

[treeProcessor:overlap]
type = OverlapRemover
listPath = eventLists
EEProcessors = lowPtSelector
EMuProcessors = lowPtSelector
MuMuProcessors = lowPtSelector
ETauProcessors = lowPtSelector
MuTauProcessors = lowPtSelector
TauTauProcessors = lowPtSelector

    """
    
    def setUp(self):
        #from ConfigParser import ConfigParser
        from helpers import BetterConfigParser as ConfigParser
        configFile = open("treePostprocessor.unittest.General.ini","w")
        configFile.write(self.configString)
        configFile.close()
        
        self.config = ConfigParser()
        self.config.read("treePostprocessor.unittest.General.ini")

    
    def tearDown(self):
        from os import remove
        try: remove("treePostprocessor.unittest.General.ini")
        except: pass    
    
    def testProducer(self):
        basePath = self.config.get("general","basePath")
        producers = getProducers(self.config, [basePath] )
        for p in producers:
            #print p.name, p.inputPaths
            if not p.name == "MergedData":                
                p.produce()

        #main(["unittest","-C","treePostprocessor.unittest.General.ini"])
        
        

