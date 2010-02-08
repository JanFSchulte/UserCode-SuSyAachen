import FWCore.ParameterSet.Config as cms

triggerResultsCounter = cms.EDAnalyzer('TriggerResultsCounter',
  prefix = cms.string("filterPathFor"),
  triggerTag = cms.InputTag("TriggerResults"),
  count = cms.VPSet(
        cms.PSet(
            name = cms.string("Same-Sign Di-Lepton Events"),
            triggerNames = cms.vstring("SSEE","SSEMu","SSMuMu","SSETau","SSMuTau","SSTauTau")
            )
        )
)

def makeFilterPaths(process, prefix="filterPathFor"):
    if process.schedule_() == None:
        process.schedule = cms.Schedule()
        for pathName in process.paths:
            process.schedule.append( process.__getattribute__(pathName) )
            
    for filter in [ i for i in process.schedule_().moduleNames() if(i in process.filters)]:
        process.__setattr__("%s%s"%(prefix,filter),cms.Path(process.__getattribute__(filter)))
        process.schedule_().append(process.__getattribute__("filterPathFor%s"%filter))
	
    for endpath in process.endpaths:
        process.schedule_().append( process.endpaths[endpath])
