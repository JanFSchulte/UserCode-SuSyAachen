import FWCore.ParameterSet.Config as cms


eeTriggerNames2017=cms.untracked.vstring( 
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_DoubleEle33_CaloIdL_MW_v",
    "HLT_DoubleEle25_CaloIdL_MW_v",
)
   
emTriggerNames2017=cms.untracked.vstring(
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",       
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",       
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",            
    "HLT_Mu27_Ele37_CaloIdL_MW_v",       
    "HLT_Mu37_Ele27_CaloIdL_MW_v",  
)
   
mmTriggerNames2017=cms.untracked.vstring(                
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", # from Run2017C
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", # only Run2017B
    "HLT_Mu37_TkMu27_v",        
)
   
htTriggerNames2017=cms.untracked.vstring(
    "HLT_PFHT180_v",
    "HLT_PFHT250_v",
    "HLT_PFHT370_v",
    "HLT_PFHT430_v",
    "HLT_PFHT510_v", 
    "HLT_PFHT590_v",
    "HLT_PFHT680_v",
    "HLT_PFHT780_v",
    "HLT_PFHT890_v",
    "HLT_PFHT1050_v" 
)     

metTriggerNames2017=cms.untracked.vstring(
    "HLT_PFMET120_PFMHT120_IDTight_v",
    "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",
)


eeTriggerNames2018=cms.untracked.vstring( 
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_DoubleEle25_CaloIdL_MW_v",
)
   
emTriggerNames2018=cms.untracked.vstring(
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",          
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",            
    "HLT_Mu27_Ele37_CaloIdL_MW_v",       
    "HLT_Mu37_Ele27_CaloIdL_MW_v",  
)
   
mmTriggerNames2018=cms.untracked.vstring(                
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", 
    "HLT_Mu37_TkMu27_v",        
)
   
htTriggerNames2018=cms.untracked.vstring(
    "HLT_PFHT180_v",
    "HLT_PFHT250_v",
    "HLT_PFHT370_v",
    "HLT_PFHT430_v",
    "HLT_PFHT510_v", 
    "HLT_PFHT590_v",
    "HLT_PFHT680_v",
    "HLT_PFHT780_v",
    "HLT_PFHT890_v",
    "HLT_PFHT1050_v" 
)     

metTriggerNames2018=cms.untracked.vstring(
    "HLT_PFMET120_PFMHT120_IDTight_v",
    "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",
)

