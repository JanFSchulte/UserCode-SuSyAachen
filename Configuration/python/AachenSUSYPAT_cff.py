import FWCore.ParameterSet.Config as cms

def reduceEventsize(process):
	#---- taus
    tauReduction(process, process.patTausPF)
    process.cleanPatTausPF.preselection = "pt > 10"
    tauReduction(process, process.patTaus)
    process.cleanPatTaus.preselection = 'pt > 10'
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
    setattr(tauCollection, "embedLeadTrack", False )
    setattr(tauCollection, "embedLeadPFCand", False )
    setattr(tauCollection, "embedSignalPFChargedHadrCands", False)
    setattr(tauCollection, "embedIsolationPFGammaCands", False)
    setattr(tauCollection, "embedSignalPFGammaCands", False)
    setattr(tauCollection, "embedIsolationPFCands", False)
    setattr(tauCollection, "embedSignalPFCands", False)
    setattr(tauCollection, "embedSignalTracks", False)
    setattr(tauCollection, "embedIsolationPFNeutralHadrCands", False)
    setattr(tauCollection, "embedIsolationPFChargedHadrCands", False)
    setattr(tauCollection, "embedIsolationTracks", False)
    setattr(tauCollection, "embedSignalPFNeutralHadrCands", False)
    setattr(tauCollection, "embedLeadPFChargedHadrCand", False)
    setattr(tauCollection, "embedLeadPFNeutralCand", False)
