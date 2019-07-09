import FWCore.ParameterSet.Config as cms

eeTriggerNames2016 = cms.untracked.vstring( 
              "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",  
              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
              "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",   
              "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v"
)
mmTriggerNames2016 = cms.untracked.vstring(  
              "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",           
              "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",            
              "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",           
              "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",            
              "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",          
              "HLT_Mu27_TkMu8_v",                          
              "HLT_Mu30_TkMu11_v"
)
emTriggerNames2016 = cms.untracked.vstring(
              "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",  
              "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",   
              "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v",   
              "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",  
              "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",  
              "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",            
              "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",         
              "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",            
              "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",        
              "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",        
              "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v",          
              "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v"
)
htTriggerNames2016 = cms.untracked.vstring(
              "HLT_PFHT125_v",
              "HLT_PFHT200_v",
              "HLT_PFHT250_v",
              "HLT_PFHT300_v",
              "HLT_PFHT350_v", 
              "HLT_PFHT400_v",
              "HLT_PFHT475_v",
              "HLT_PFHT600_v",
              "HLT_PFHT650_v",
              "HLT_PFHT800_v",
              "HLT_PFHT900_v",
              
)

metTriggerNames2016=cms.untracked.vstring(
    "HLT_PFMET120_PFMHT90_IDTight_v",
    "HLT_PFMET120_PFMHT100_IDTight_v",
    "HLT_PFMET120_PFMHT110_IDTight_v",
    "HLT_PFMET120_PFMHT120_IDTight_v",
)

eeTriggerNames2017=cms.untracked.vstring( 
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_DoubleEle33_CaloIdL_MW_v",
    "HLT_DoubleEle25_CaloIdL_MW_v",
    "HLT_DoubleEle27_CaloIdL_MW_Edge",
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

