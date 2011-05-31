#!/usr/bin/env python
'''
Created on 26.05.2011

@author: heron
'''

class TreeProcessor:
    def __init__(self, config, name):
        self.config = config
        self.section = "treeProcessor:%s"%name
        self.nThEntry = 0
        self.nEntries = None
        
    def prepareSrc(self, tree, object):
        self.nEntries = tree.GetEntries()
        
    
    def prepareDest(self, tree, object):
        pass

    def processEvent(self, event, object):
        import sys
        self.nThEntry +=1
        if self.nThEntry%(self.nEntries*0.001) < 1: 
            sys.stdout.write("\r%.1f%%" %(self.nThEntry*1./(self.nEntries*0.01)))    
            sys.stdout.flush()  
        return True
    
class SimpleSelector(TreeProcessor):
    def __init__(self, config, name):
        TreeProcessor.__init__(self, config, name)
        
    def processEvent(self, event, object):
        return event.ht > 300
    
class SimpleWeighter(TreeProcessor):
    def __init__(self, config, name):
        TreeProcessor.__init__(self, config, name)
        self.weight = {}
        
    def prepareSrc(self, src, object):
        #src.SetBranchStatus("weight", 0)
        pass
        
    def prepareDest(self, dest, object):
        from array import array
        self.weight[object] = array("f",[1.0])
        dest.Branch("weight2",self.weight[object],"weight2/F")

    def processEvent(self, event, object):
        self.weight[object][0] = event.weight*(event.pt1+event.pt2)
        return True

class FakeWeighter(TreeProcessor):
    def __init__(self, config, name):
        from plotFakeRate import parsePSet
        TreeProcessor.__init__(self, config, name)
        frPath = config.get(self.section, "fakePSet")
        self.pSet = parsePSet(frPath)
        self.ptMax = 74.9
        self.weight = {}
        
    def prepareSrc(self, src, object):
        TreeProcessor.prepareSrc(self, src, object)
        #src.SetBranchStatus("weight", 0)
        pass
        
    def prepareDest(self, dest, object):
        from array import array
        TreeProcessor.prepareDest(self, dest, object)
        self.weight[object] = array("f",[-1.0])
        dest.Branch("fakeWeight",self.weight[object],"fakeWeight/F")
        
    def processEvent(self, event, object):
        TreeProcessor.processEvent(self, event, object)
        self.weight[object][0] = -1
        if event.tauDiscr < 0.5:
            pt = event.pt
            if pt > self.ptMax: pt =  self.ptMax
            f = self.pSet["center"].fakeRate({"pt":pt, "eta":event.eta})
            self.weight[object][0] = f/(1-f)
        return True

    
class OverlapRemover(TreeProcessor):
    def __init__(self, config, name):
        TreeProcessor.__init__(self, config, name)
        self.keepEvents = {}
    
    def prepareSrc(self, src, object):
        for ev in src:
            fingerPrint = (ev.runNr, ev.lumiSec, ev.eventNr) 
            if not fingerPrint in self.keepEvents or self.keepEvents[fingerPrint][1]< ev.pt1+ev.pt2:
                self.keepEvents[fingerPrint] = (object, ev.pt1+ev.pt2) 
    
    def processEvent(self, event, object):
        fingerPrint = (event.runNr, event.lumiSec, event.eventNr)        
        #if fingerPrint in self.keepEvents and not object == self.keepEvents[fingerPrint][0]:
        #    print "skipping", fingerPrint, object, event.pt1+event.pt2,   self.keepEvents[fingerPrint]
        return fingerPrint in self.keepEvents and object == self.keepEvents[fingerPrint][0] 

class TreeProducer:
    def __init__(self, config, inputPaths, name = None):
        from os.path import split as splitPath
        from os.path import exists as pathExists
        from os import makedirs
        from sys import modules
        if name == None:
            assert len(inputPaths) == 1, "can only determine names for singlePaths automatically. Got '%s'"%inputPaths
            name = splitPath(inputPaths[0])[1].split(".")[2]
        self.config = config
        self.name = name
        self.tasks = list(set([splitPath(i)[1].split(".")[1] for i in inputPaths]))
        self.flags = list(set([splitPath(i)[1].split(".")[0] for i in inputPaths]))
        self.inputPaths = inputPaths
        self.outPath = config.get("general","outPath")
        if not pathExists(self.outPath):
            makedirs(self.outPath)
        self.treeProcessors = {}
        for section in self.config.sections():
            if section.startswith("treeProcessor:"):
                processorName = section.split("treeProcessor:")[1]
                processorType = self.config.get(section,"type")
                #warning black magic ahead :P
                self.treeProcessors[processorName] = getattr(modules[globals()["__name__"]],processorType)(self.config, processorName)
        
    def produce(self):
        from ROOT import TFile, TChain
        outFile = TFile("%s/%s.%s.%s.root"%(self.outPath, "".join(self.flags), "processed" , self.name),"recreate")
        for section in self.config.sections():
            trees = None
            if section.startswith("dileptonTree:"):
                trees = self._getDileptonTrees(section)
                treeName = "DileptonTree"
                subDirName = section.split("dileptonTree:")[1]
            if section.startswith("isoTree:"):
                trees = self._getIsoTrees(section)
                treeName = "Iso"
                subDirName = section.split("isoTree:")[1]
            if not trees == None:
                outDir = None 
                for object in trees:
                    srcTree = TChain("%s%s"%(object, treeName))
                    processors = self.config.get(section,"%sProcessors"%object).split()
                    filter = " and ".join(processors)
                    if self.config.has_option(section,"%sFilter"%object):
                        filter = self.config.get(section,"%sFilter"%object)
                    for treePath in trees[object]:
                        #srcFile = TFile(filePath,"r")
                        #srcTree = srcFile.Get(treePath)
                        srcTree.Add(treePath)
                        print "adding", treePath
                    srcTree.SetBranchStatus("*", 1)
                    for processorName in processors:
                        self.treeProcessors[processorName].prepareSrc(srcTree, object)
                    if not outDir:
                        outDir = outFile.mkdir(subDirName)
                    outFile.cd(subDirName)
                    destTree = srcTree.CloneTree(0)
                    for processorName in processors:
                        self.treeProcessors[processorName].prepareDest(destTree, object)
                    for i in srcTree:
                        processingResults = {}
                        for processorName in processors:
                            processingResults[processorName] = self.treeProcessors[processorName].processEvent(srcTree, object)
                        if eval(filter, processingResults):
                            destTree.Fill()
                    #srcFile.Close()
                    outFile.Write()
                #from pprint import pprint
                #pprint( trees)
        
        outFile.Close()
                
                
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
            datasetPaths = []
            datasetExpression = self.config.get(section, "%sDataset"%object)
            datasetSelection = self.config.get(section, "%sSelection"%object)
            for path in self.inputPaths:
                if not match(datasetExpression, splitPath(path)[1].split(".")[2]) == None:
                    datasetPaths.append(path)
            if datasetPaths == []:
                datasetSelection = self.config.get(section, "otherSelection")
                datasetPaths = self.inputPaths
            
            result[object] = []
            for path in datasetPaths:
                task = None
                for t in self.tasks:
                    if ".%s."%t in splitPath(path)[1]:
                        assert task == None, "unable to disambiguate tasks '%s' matches both '%s' and '%s'"(path, task, t)
                        task = t
                result[object].append( "%s/%s%s%s/%sDileptonTree"%(path,task,datasetSelection,treeProducerName,object))
        return result
        
def getProducers(config, path):
    from glob import glob
    from re import match
    from os.path import split as splitPath
    mcExpression = config.get("general","MCDatasets")
    result = []
    dataPaths = []
    for inputPath in glob("%s/*.root"%path):
        if match(mcExpression, splitPath(inputPath)[1]) ==None:
            dataPaths.append(inputPath)
        else:   
            result.append(TreeProducer(config, [inputPath]))
    if len(dataPaths) > 0:
        result.append(TreeProducer(config, dataPaths, "MergedData"))
    return result

    
    
def main(argv = None):
    import sys
    from ConfigParser import ConfigParser
    from optparse import OptionParser
    if argv == None:
        argv = sys.argv[1:]
    parser = OptionParser()
    parser.add_option("-C", "--Config", dest="Config", action="append", default=[],
                          help="Main configuration file. Can be given multiple times in case of split configurations. Default is Input/default.ini")
    (opts, args) = parser.parse_args(argv)

    if opts.Config == []:
        opts.Config = [ "default.ini" ]
    
    config = ConfigParser()
    config.read(opts.Config)
    
    basePath = config.get("general","basePath")
    producers = getProducers(config, basePath)
    for p in producers:
        p.produce()
    
if __name__ == '__main__':
    main()

import unittest    

class plotTest(unittest.TestCase):
    configString ="""
[general]
tasks = pfBaseCuts pfBaseCutsSingleMu pfBaseCutsSingleE pfBaseCutsTauPlusX
basePath = /Users/heron/Documents/superSymmetry/results/diLeptonTaus/diLeptons41x/susy0447v5/input
MCDatasets = .*_Spring11
outPath = processedTrees

[dileptonTree NoCuts]
treeProducerName = TaNCTrees
objects = EE EMu MuMu ETau MuTau TauTau
EEDataset = SingleElectron.* 
EESelection = HLTE
EEProcessors = htSelector ptSumWeighter overlap
EEFilter = True
EMuDataset = SingleMu_.*
EMuSelection = HLTIsoMu
EMuProcessors = htSelector ptSumWeighter overlap
EMuFilter = htSelector 
MuMuDataset = SingleMu_.*
MuMuSelection = HLTMu
MuMuProcessors = htSelector overlap
ETauDataset = TauPlusX_.*
ETauSelection = HLTETau
ETauProcessors = htSelector overlap
MuTauDataset = TauPlusX_.*
MuTauSelection = HLTMuTau
MuTauProcessors = htSelector overlap
TauTauDataset = TauPlusX_.*
TauTauSelection = HLTTauTauHT
TauTauProcessors = htSelector overlap
OtherSelection = 

[isoTree:NoCuts]
treeProducerName = TnPTaNCTauTrees
Dataset = TauPlusX
Selection = 
Processors = tauFakeWeights
PtherSelection = 

[treeProcessor:htSelector]
type = SimpleSelector

[treeProcessor:ptSumWeighter]
type = SimpleWeighter

[treeProcessor:tauFakeWeights]
type = FakeWeighter
fakePSet = /Users/heron/Documents/superSymmetry/results/diLeptonTaus/diLeptons41x/susy0447v5/tauDataHT_cff.py
selection = tauDiscr > 0.5

[treeProcessor:overlap]
type = OverlapRemover
    """
    
    def setUp(self):
        from ConfigParser import ConfigParser
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
        producers = getProducers(self.config, basePath)
        for p in producers:
            #print p.name, p.inputPaths
            if p.name == "MergedData":                
                p.produce()

        #main(["unittest","-C","treePostprocessor.unittest.General.ini"])
        
        
