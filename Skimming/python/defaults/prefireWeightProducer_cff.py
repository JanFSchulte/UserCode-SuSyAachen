import FWCore.ParameterSet.Config as cms

def prefireWeightProducer(process):
        process.prefiringweight = cms.EDProducer("L1ECALPrefiringWeightProducer",
                                 ThePhotons = cms.InputTag("slimmedPhotons"),
                                 TheJets = cms.InputTag("slimmedJets"),
                                 L1Maps = cms.string("${CMSSW_BASE}/src/SuSyAachen/DiLeptonHistograms/data/L1PrefiringMaps_new.root"), # update this line with the location of this file
                                 DataEra = cms.string("2017BtoF"), #Use 2016BtoH for 2016
                                 UseJetEMPt = cms.bool(False), #can be set to true to use jet prefiring maps parametrized vs pt(em) instead of pt
                                 PrefiringRateSystematicUncty = cms.double(0.2) #Minimum relative prefiring uncty per object
                                 )

        
        process.seqprefireWeightProducer = cms.Sequence(process.prefiringweight)
        process.prefireWeightPath = cms.Path(process.seqprefireWeightProducer)
