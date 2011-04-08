import FWCore.ParameterSet.Config as cms

def reduceEventsize(process):
	#---- taus
    tauReduction(process, process.patTausPF)
    process.pfTausPF.discriminators = cms.VPSet()
    process.cleanPatTausPF.preselection = "1"
    process.selectedPatTausPF.cut = "pt > 10 & tauID('leadingTrackFinding') > 0.5"
    
    tauReduction(process, process.patTaus)
    process.pfTaus.discriminators = cms.VPSet()
    process.cleanPatTaus.preselection = "pt > 10 & tauID('leadingTrackFinding') > 0.5"
    
    if hasattr(process,"pfTausPFShrinkingCone"):
	    tauReduction(process, process.patTausPFShrinkingCone)
	    process.pfTausPFShrinkingCone.discriminators = cms.VPSet()
	    process.selectedPatTausPFShrinkingCone.cut = "pt > 10 & tauID('leadingTrackFinding') > 0.5"
    	
    
    # for later
    #    process.cleanPatTaus.preselection = 'pt > 10 & tauID("againstMuon") > 0.5 & tauID("againstElectron") > 0.5 & tauID("byTaNCfrOnePercent") > 0.5'
    #---- jets
    process.patJetsAK5PF.embedGenJetMatch = True
    process.patJetsAK5PF.embedPFCandidates = False
    process.patJetsPF.embedGenJetMatch = True
    process.patJetsPF.embedPFCandidates = False
    process.patJets.embedGenJetMatch = True
    process.patJets.embedPFCandidates = False
    process.patJets.embedCaloTowers = False
    


def tauReduction(process, tauCollection):    
    setattr(tauCollection, "isoDeposits", cms.PSet() )
    setattr(tauCollection, "embedLeadTrack", True )
    setattr(tauCollection, "embedLeadPFCand", False )
    setattr(tauCollection, "embedSignalPFChargedHadrCands", False)
    setattr(tauCollection, "embedIsolationPFGammaCands", False)
    setattr(tauCollection, "embedSignalPFGammaCands", False)
    setattr(tauCollection, "embedIsolationPFCands", False)
    setattr(tauCollection, "embedSignalPFCands", False)
    setattr(tauCollection, "embedSignalTracks", True)
    setattr(tauCollection, "embedIsolationPFNeutralHadrCands", False)
    setattr(tauCollection, "embedIsolationPFChargedHadrCands", False)
    setattr(tauCollection, "embedIsolationTracks", False)
    setattr(tauCollection, "embedSignalPFNeutralHadrCands", False)
    setattr(tauCollection, "embedLeadPFChargedHadrCand", False)
    setattr(tauCollection, "embedLeadPFNeutralCand", False)
    
def additionalTaus(process,postfix="PF"):
    from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet, massSearchReplaceAnyInputTag
    cloneProcessingSnippet(process,getattr(process,"pfTausBaseSequence"+postfix),"ShrinkingCone")
    cloneProcessingSnippet(process,getattr(process,"makePatTaus"+postfix),"ShrinkingCone")
    cloneProcessingSnippet(process,getattr(process,"patShrinkingConePFTauDiscrimination"+postfix),"ShrinkingCone")
    
    setattr(process,"selectedPatTaus"+postfix+"ShrinkingCone",
            getattr(process,"selectedPatTaus"+postfix).clone())

    setattr(process,"pfTaus"+postfix+"ShrinkingCone",
            getattr(process,"pfTaus"+postfix).clone())
    

    massSearchReplaceAnyInputTag(getattr(process,"patShrinkingConePFTauDiscrimination"+postfix+"ShrinkingCone"), cms.InputTag("pfTaus"+postfix), cms.InputTag("pfTaus"+postfix+"ShrinkingCone"))
    massSearchReplaceAnyInputTag(getattr(process,"patPFTauIsolation"+postfix+"ShrinkingCone"), cms.InputTag("pfTaus"+postfix), cms.InputTag("pfTaus"+postfix+"ShrinkingCone"))
    
#    getattr(process,"patTaus"+postfix+"ShrinkingCone").addDecayMode = True
    massSearchReplaceAnyInputTag(getattr(process,"makePatTaus"+postfix+"ShrinkingCone"),cms.InputTag("pfTaus"+postfix), cms.InputTag("pfTaus"+postfix+"ShrinkingCone"))
#    getattr(process,"patTaus"+postfix+"ShrinkingCone").decayModeSrc = "shrinkingConePFTauDecayModeProducer"+postfix+"ShrinkingCone"
    getattr(process,"selectedPatTaus"+postfix+"ShrinkingCone").src = cms.InputTag("patTaus"+postfix+"ShrinkingCone")
    getattr(process,"pfTaus"+postfix+"ShrinkingCone").src = cms.InputTag("pfTausBase"+postfix+"ShrinkingCone")
    
    from PhysicsTools.PatAlgos.tools.pfTools import adaptPFTaus
    adaptPFTaus(process,"hpsPFTau",postfix=postfix)
    
    getattr(process,"pfTauSequence"+postfix).replace(
        getattr(process,"pfTausBaseSequence"+postfix),
        getattr(process,"pfTausBaseSequence"+postfix) +
        getattr(process,"pfTausBaseSequence"+postfix+"ShrinkingCone") 
        )

    getattr(process,"pfTauSequence"+postfix).replace(
        getattr(process,"pfTaus"+postfix),
        getattr(process,"pfTaus"+postfix) +
        getattr(process,"pfTaus"+postfix+"ShrinkingCone") 
        )

	#not used anymore
    #getattr(process,"patCandidates"+postfix).replace(
    #    getattr(process,"makePatTaus"+postfix),
    #    getattr(process,"makePatTaus"+postfix) +
    #    getattr(process,"makePatTaus"+postfix+"ShrinkingCone") 
    #    )

    getattr(process,"patDefaultSequence"+postfix).replace(
        getattr(process,"patShrinkingConePFTauDiscrimination"+postfix),
        getattr(process,"patShrinkingConePFTauDiscrimination"+postfix) +
        getattr(process,"patShrinkingConePFTauDiscrimination"+postfix+"ShrinkingCone") 
        )

    for modName in getattr(process,"makePatTaus"+postfix).moduleNames():
        if modName in getattr(process,"patDefaultSequence"+postfix).moduleNames():
            getattr(process,"patDefaultSequence"+postfix).replace(
                getattr(process,modName),
                getattr(process,modName) +
                getattr(process,modName+"ShrinkingCone") 
    	    	)

    getattr(process,"patDefaultSequence"+postfix).replace(
        getattr(process,"selectedPatTaus"+postfix),
        getattr(process,"selectedPatTaus"+postfix) +
        getattr(process,"selectedPatTaus"+postfix+"ShrinkingCone") 
        )

    #Relax tau preselection
#    getattr(process,"selectedPatTaus"+postfix).cut = 'pt > 15 & abs(eta) < 2.5'
    getattr(process,"pfTaus"+postfix).discriminators = cms.VPSet()
    
    getattr(process,"selectedPatTaus"+postfix+"ShrinkingCone").cut = 'pt > 10'
    getattr(process,"pfTaus"+postfix+"ShrinkingCone").discriminators = cms.VPSet()
	#discriminators
    
    setattr(process,"hpsPFTauDiscriminationAgainstElectron2D"+postfix,
            getattr(process,"hpsPFTauDiscriminationAgainstElectron"+postfix).clone(
        ApplyCut_ElectronPreID_2D = cms.bool(True),
        ApplyCut_PFElectronMVA =  cms.bool(False)
        )
            )
    setattr(process,"hpsPFTauDiscriminationAgainstElectronCrackRem"+postfix,
            getattr(process,"hpsPFTauDiscriminationAgainstElectron"+postfix).clone(
        ApplyCut_EcalCrackCut = cms.bool(True),
        ApplyCut_PFElectronMVA =  cms.bool(False)
        )
            )
    
    setattr(process,"shrinkingConePFTauDiscriminationAgainstElectron2D"+postfix+"ShrinkingCone",
            getattr(process,"shrinkingConePFTauDiscriminationAgainstElectron"+postfix+"ShrinkingCone").clone(
        ApplyCut_ElectronPreID_2D = cms.bool(True),
        ApplyCut_PFElectronMVA =  cms.bool(False)
        )
            )
    setattr(process,"shrinkingConePFTauDiscriminationAgainstElectronCrackRem"+postfix+"ShrinkingCone",
            getattr(process,"shrinkingConePFTauDiscriminationAgainstElectron"+postfix+"ShrinkingCone").clone(
        ApplyCut_EcalCrackCut = cms.bool(True),
        ApplyCut_PFElectronMVA =  cms.bool(False)
        )
            )
	
    s = getattr(process,"patHPSPFTauDiscrimination"+postfix) 
    s +=     getattr(process,"hpsPFTauDiscriminationAgainstElectron2D"+postfix)
    s = getattr(process,"patHPSPFTauDiscrimination"+postfix) 
    s +=     getattr(process,"hpsPFTauDiscriminationAgainstElectronCrackRem"+postfix)

    s= getattr(process,"patShrinkingConePFTauDiscrimination"+postfix+"ShrinkingCone") 
    s +=     getattr(process,"shrinkingConePFTauDiscriminationAgainstElectron2D"+postfix+"ShrinkingCone")
    s= getattr(process,"patShrinkingConePFTauDiscrimination"+postfix+"ShrinkingCone") 
    s +=     getattr(process,"shrinkingConePFTauDiscriminationAgainstElectronCrackRem"+postfix+"ShrinkingCone")

    getattr(process,"patTaus"+postfix).tauIDSources = cms.PSet(
        leadingTrackFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"+postfix),
        byLooseIsolation = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"+postfix),
        byMediumIsolation = cms.InputTag("hpsPFTauDiscriminationByMediumIsolation"+postfix),
        byTightIsolation = cms.InputTag("hpsPFTauDiscriminationByTightIsolation"+postfix),
        againstElectron = cms.InputTag("hpsPFTauDiscriminationAgainstElectron"+postfix),
        againstElectron2D = cms.InputTag("hpsPFTauDiscriminationAgainstElectron2D"+postfix),
        againstElectronCrackRem = cms.InputTag("hpsPFTauDiscriminationAgainstElectronCrackRem"+postfix),
        againstMuon = cms.InputTag("hpsPFTauDiscriminationAgainstMuon"+postfix)
        )
    getattr(process,"patTaus"+postfix+"ShrinkingCone").tauIDSources = cms.PSet(
        leadingTrackFinding = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackFinding"+postfix+"ShrinkingCone"),
        leadingTrackPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackPtCut"+postfix+"ShrinkingCone"),
        leadingPionPtCut = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingPionPtCut"+postfix+"ShrinkingCone"),
        trackIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolation"+postfix+"ShrinkingCone"),
        trackIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion"+postfix+"ShrinkingCone"),
        ecalIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolation"+postfix+"ShrinkingCone"),
        ecalIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion"+postfix+"ShrinkingCone"),
        byIsolation = cms.InputTag("shrinkingConePFTauDiscriminationByIsolation"+postfix+"ShrinkingCone"),
        byIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByIsolationUsingLeadingPion"+postfix+"ShrinkingCone"),
        againstElectron = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectron"+postfix+"ShrinkingCone"),
        againstElectron2D = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectron2D"+postfix+"ShrinkingCone"),
        againstElectronCrackRem = cms.InputTag("shrinkingConePFTauDiscriminationAgainstElectronCrackRem"+postfix+"ShrinkingCone"),
        againstMuon = cms.InputTag("shrinkingConePFTauDiscriminationAgainstMuon"+postfix+"ShrinkingCone"),
        byTaNC = cms.InputTag("shrinkingConePFTauDiscriminationByTaNC"+postfix+"ShrinkingCone"),
        byTaNCfrOnePercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrOnePercent"+postfix+"ShrinkingCone"),
        byTaNCfrHalfPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrHalfPercent"+postfix+"ShrinkingCone"),
        byTaNCfrQuarterPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrQuarterPercent"+postfix+"ShrinkingCone"),
        byTaNCfrTenthPercent = cms.InputTag("shrinkingConePFTauDiscriminationByTaNCfrTenthPercent"+postfix+"ShrinkingCone")

        )
    
