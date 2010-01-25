import FWCore.ParameterSet.Config as cms
def GenLeptons(process):
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")# import HepPDTESSource

    process.electronGenParticles = cms.EDFilter( "GenParticlePruner", 
	    src = cms.InputTag("genParticles"),
	    select = cms.vstring(
	    "drop  *", # this is the default
	    "keep (pdgId = {e+} &status=1)| (pdgId = {e-}  &status=1)",
	    )
    )

    process.muonGenParticles = cms.EDFilter( "GenParticlePruner", 
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
        "drop  *", # this is the default
        "keep (pdgId = {mu+} & status=1) | (pdgId = {mu-} & status=1)",
    #    "drop (pdgId = {mu+} & status = 2) | (pdgId = {mu-} & status = 2)"
        )
    )

    process.genTausWithHistory = cms.EDFilter( "GenParticlePruner", 
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
        "drop  *", # this is the default
        "keep++ pdgId = {tau+} | pdgId = {tau-}",
        #NOTE the taus will decay one needs to select the decay mode and then drop the daughters
        )
    )

    process.hadronicGenTaus = cms.EDFilter("GenDaughterExcluder", 
        src = cms.InputTag("genTausWithHistory"),
        daughterIds = cms.vint32( 11,-11, 13,-13) #e+ e- mu+ mu-
    )
    
    process.tauGenParticles = cms.EDFilter( "GenParticlePruner", 
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
        "drop  *", # this is the default
        "keep pdgId = {tau+} | pdgId = {tau-}",
        "drop  status = 3"
        )
    )
    
    process.hadronicTauGenParticles = cms.EDFilter( "GenParticlePruner", 
        src = cms.InputTag("hadronicGenTaus"),
        select = cms.vstring(
        "drop  *", # this is the default
        "keep pdgId = {tau+} | pdgId = {tau-}",
        "drop  status = 3"
        )
    )


########## --- Prompt selections
    from SuSyAachen.Skimming.GenPromptSelector_cfi import GenPromptSelector

    process.genEWithHistory = cms.EDFilter( "GenParticlePruner", 
	    src = cms.InputTag("genParticles"),
	    select = cms.vstring(
	    "drop  *", # this is the default
	    "++keep (pdgId = {e+} &status=1)| (pdgId = {e-}  &status=1)",
	    )
    )
    process.promptElectrons = GenPromptSelector.clone(
        src = "genEWithHistory"
    )

    process.genMuWithHistory = cms.EDFilter( "GenParticlePruner", 
	    src = cms.InputTag("genParticles"),
	    select = cms.vstring(
	    "drop  *", # this is the default
	    "++keep (pdgId = {mu+} &status=1)| (pdgId = {mu-}  &status=1)",
	    )
    )
    process.promptMuons = GenPromptSelector.clone(
        src = "genMuWithHistory"
    )

    #process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
    #    src = cms.InputTag("genEWithHistory"),#genTausWithHistory"),#tauGenParticles"),#genParticles"),#
    #    printP4 = cms.untracked.bool(False),
    #    printPtEtaPhi = cms.untracked.bool(False),
    #    printVertex = cms.untracked.bool(False),    
    #    printStatus = cms.untracked.bool(True),
    #    printIndex = cms.untracked.bool(False),
    #    status = cms.untracked.vint32( 1, 2, 3 )
    #  ) 

    process.seqGenLeptons = cms.Sequence(process.electronGenParticles + process.muonGenParticles + process.tauGenParticles
                                         + (process.genTausWithHistory * process.hadronicGenTaus * process.hadronicTauGenParticles)
                                         + (process.genEWithHistory * process.promptElectrons) + (process.genMuWithHistory * process.promptMuons)
                                         )
