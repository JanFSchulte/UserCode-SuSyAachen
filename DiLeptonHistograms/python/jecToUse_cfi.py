import FWCore.ParameterSet.Config as cms
# central place to edit used JEC so it is always consistent
# don't forget to also add/remove the DB files to/from additional input
# files for crab


# whether to load from GT or from DB file
dataUseDB = {}
dataUseDB["2016"] = True
dataUseDB["2017"] = True
dataUseDB["2018"] = True

# db name, can be empty string if not used
dataDB = {}
dataDB["2016"] = "Summer16_07Aug2017All_V11_DATA"
dataDB["2017"] = "Fall17_17Nov2017_V32_94X_DATA"
dataDB["2018"] = "Autumn18_RunABCD_V19_DATA"

# whether to load from GT or from DB file
mcUseDB = {}
mcUseDB["2016"] = True
mcUseDB["2017"] = True
mcUseDB["2018"] = True

# db name, can be empty string if not used
mcDB = {}
mcDB["2016"] = "Summer16_07Aug2017_V11_MC"
mcDB["2017"] = "Fall17_17Nov2017_V32_94X_MC"
mcDB["2018"] = "Autumn18_V19_MC"

fastSimUseDB = {}
fastSimUseDB["2016"] = True
fastSimUseDB["2017"] = True
fastSimUseDB["2018"] = True

fastSimDB =  {}
fastSimDB["2016"] = "Spring16_25nsFastSimV1_MC"
fastSimDB["2017"] = "Fall17_FastSimV1_MC"
fastSimDB["2018"] = "Autumn18_FastSimV1_MC"
