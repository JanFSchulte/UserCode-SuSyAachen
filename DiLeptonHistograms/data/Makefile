#doku in https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
ANA_JSON=/home/home4/institut_1b/edelhoff/superSymmetry/SubmitScripts/Input/Cert_190456-195016_8TeV_PromptReco_Collisions12_JSON.txt
#/home/home4/institut_1b/edelhoff/superSymmetry/SubmitScripts/Input/Cert_190456-193557_8TeV_PromptReco_Collisions12_JSON.txt
#/home/home4/institut_1b/edelhoff/superSymmetry/SubmitScripts/Input/Cert_190456-191276_8TeV_PromptReco_Collisions12_JSON.txt
PU_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-195165.txt
#5fbi? /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-196531_patch2.txt
#1p6fbi/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-193575.txt
#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-191276.txt
PU_CALC=pileupCalc.py --calcMode true --minBiasXsec 69300 --maxPileupBin 60 --numPileupBins 60 

mc_PU_3D_dist.root:
	ls

mc_PU_dist.root:
	ls

PU_dist.root:
	$(PU_CALC) -i $(ANA_JSON) --inputLumiJSON $(PU_JSON) --calcMode observed $@

PU_3D_dist.root:
	ls
