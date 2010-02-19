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

    process.genTausWithFuture = cms.EDFilter( "GenParticlePruner", 
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
        "drop  *", # this is the default
        "keep++ pdgId = {tau+} | pdgId = {tau-}",
        #NOTE the taus will decay one needs to select the decay mode and then drop the daughters
        )
    )

    process.hadronicGenTaus = cms.EDFilter("GenDaughterExcluder", 
        src = cms.InputTag("genTausWithFuture"),
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

    process.promptParticles = GenPromptSelector.clone(
        src = "genParticles"
        )
 
    process.promptElectrons = cms.EDFilter( "GenParticlePruner", 
            src = cms.InputTag("promptParticles"),
            select = cms.vstring(
            "drop  *", # this is the default
            "keep pdgId = {e+} | pdgId = {e-} "
            )
    )
    process.promptMuons = cms.EDFilter( "GenParticlePruner", 
            src = cms.InputTag("promptParticles"),
            select = cms.vstring(
            "drop  *", # this is the default
            "keep pdgId = {mu+} | pdgId = {mu-}"
            )
    )
    process.promptTaus = cms.EDFilter( "GenParticlePruner", 
            src = cms.InputTag("promptParticles"),
            select = cms.vstring(
            "drop  *", # this is the default
            "keep pdgId = {tau+} | pdgId = {tau-}"
            )
    )

    process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
        src = cms.InputTag("genParticles"),#"genEWithHistory"),#genTausWithHistory"),#tauGenParticles"),#genParticles"),#
        printP4 = cms.untracked.bool(False),
        printPtEtaPhi = cms.untracked.bool(False),
        printVertex = cms.untracked.bool(False),    
        printStatus = cms.untracked.bool(True),
        printIndex = cms.untracked.bool(False),
        status = cms.untracked.vint32( 1,2, 3 )
      ) 
    process.dump = cms.EDAnalyzer('EventContentAnalyzer') 



    process.seqGenLeptons = cms.Sequence(process.electronGenParticles + process.muonGenParticles + process.tauGenParticles
                                         + (process.genTausWithFuture * process.hadronicGenTaus * process.hadronicTauGenParticles)
                                         + process.promptParticles * ( process.promptElectrons + process.promptMuons + process.promptTaus)
                                         + process.printTree
                                         )
