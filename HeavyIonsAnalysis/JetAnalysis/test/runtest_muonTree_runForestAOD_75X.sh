#!/bin/bash

##CMSSW_7_5_8, cmsRun test, 100 events for pp MC/DATA and 1000 events for PbPb DATA
##test muonTree plugin in forest scripts, where it will potentially be used

cmsRun runForestAOD_PbPb_DATA_75X.py  >  testMuonTree_AOD_PbPb_DATA_75X.log     
cmsRun runForestAOD_pp_MC_75X.py      >  testMuonTree_AOD_pp_MC_75X.log
cmsRun runForestAOD_pp_DATA_75X.py    >  testMuonTree_AOD_pp_DATA_75X.log

##no good test file @ RU T3, leaving out for now
#cmsRun runForestAOD_PbPb_MB_75X.py      testMuonTree_AOD_PbPb_MB_75X.log
#cmsRun runForestAOD_PbPb_MIX_75X.py     testMuonTree_AOD_PbPb_MIX_75X.log
