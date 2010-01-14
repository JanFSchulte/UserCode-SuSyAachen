import FWCore.ParameterSet.Config as cms

def SuperClusters(process):
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    #  SuperClusters  ################
    process.HybridSuperClusters = cms.EDProducer("ConcreteEcalCandidateProducer",
            src = cms.InputTag("correctedHybridSuperClusters"),
            particleType = cms.string('gamma')
        )
    process.EBSuperClusters = cms.EDFilter("CandViewSelector",
            src = cms.InputTag("HybridSuperClusters"),
            cut = cms.string('abs( eta ) < 1.4442')
        )

    process.EndcapSuperClusters = cms.EDProducer("ConcreteEcalCandidateProducer",
            src = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
            particleType = cms.string('gamma')
        )
    process.EESuperClusters = cms.EDFilter("CandViewSelector",
            src = cms.InputTag("EndcapSuperClusters"),
            cut = cms.string('abs( eta ) > 1.560 & abs( eta ) < 3.0')
        )

    process.allSuperClusters = cms.EDFilter("CandMerger",
            src = cms.VInputTag(cms.InputTag("EBSuperClusters"), cms.InputTag("EESuperClusters"))
        )


    process.seqSuperClusters = cms.Sequence((process.HybridSuperClusters*process.EBSuperClusters+process.EndcapSuperClusters*process.EESuperClusters)*process.allSuperClusters) 
