import FWCore.ParameterSet.Config as cms
def SUSYPARS(process):  
    process.load('SUSYBSMAnalysis.SusyCAF.SusyCAF_Scan_cfi')

    from SuSyAachen.Histograms.scanCounter_cfi import scanCounter
    from SuSyAachen.DiLeptonHistograms.DiLeptonTrees_cfi import mSugraVars
    process.scanCounter = scanCounter.clone(
        susyVars = cms.VPSet(mSugraVars),
        )
    process.seqSUSYPARS = cms.Sequence(process.susycafscan * process.scanCounter)

                         
