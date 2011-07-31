import FWCore.ParameterSet.Config as cms
def SUSYPARS(process):  
    process.load('SUSYBSMAnalysis.SusyCAF.SusyCAF_Scan_cfi')
    process.seqSUSYPARS = process.susycafscan
                         
