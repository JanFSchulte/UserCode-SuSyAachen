#
#  SUSY-PAT configuration file
#
#  PAT configuration for the SUSY group - 3X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV6
#

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#-- Meta data to be logged in DBS ---------------------------------------------
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.17 $'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/PhysicsTools/Configuration/test/SUSY_pattuple_cfg.py,v $'),
    annotation = cms.untracked.string('SUSY pattuple definition')
)

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#-- Input Source --------------------------------------------------------------
process.source.fileNames = [
#    'rfio://?svcclass=cmscafuser&path=/castor/cern.ch/cms/store/caf/user/fronga/V6production/PYTHIA6_SUSY_LM0_sftsht_10TeV_cff_py_RAW2DIGI_RECO_1.root'
'file:/user/edelhoff/mcData/CMSSW_314/QCD_Pt470-Summer09-MC_31X_V3-v1-GEN-SIM-RECO_1.root'
    ]
process.maxEvents.input = 100
# Due to problem in production of LM samples: same event number appears multiple times
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

#-- Calibration tag -----------------------------------------------------------
# Should match input file's tag
process.GlobalTag.globaltag = 'STARTUP31X_V1::All'


############################# START SUSYPAT specifics ####################################
process.load("SuSyAachen.Skimming.SUSY_pattuple_cff")

from SuSyAachen.Skimming.SUSY_pattuple_cff import addSUSYCollections
addSUSYCollections( process )

from SuSyAachen.Skimming.SUSY_pattuple_cff import getSUSY_pattuple_outputCommands
SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process.name_() )
############################## END SUSYPAT specifics ####################################


#-- Output module configuration -----------------------------------------------
process.out.fileName = 'SUSYPAT.root'       # <-- CHANGE THIS TO SUIT YOUR NEEDS

# Custom settings
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = cms.untracked.vstring('drop *', *SUSY_pattuple_outputCommands )


#-- Execution path ------------------------------------------------------------
# Full path
process.p = cms.Path( process.seqSUSY_pattuple )

