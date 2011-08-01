import FWCore.ParameterSet.Config as cms
def DumpEventContent(process):  
    process.dump = cms.EDAnalyzer('EventContentAnalyzer')
    process.seqDumpEventContent = cms.Sequence(process.dump)
                         
