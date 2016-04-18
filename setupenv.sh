##SCRAM_ARCH for specific CMSSW_XYZ area
export SCRAM_ARCH="slc6_amd64_gcc491" 

echo "sourcing condor+OSG software"
source /condor/HTCondor/current/condor.sh
source /osg/osg3.2/osg-client/setup.sh

echo 'sourcing cms software'
source "$VO_CMS_SW_DIR/cmsset_default.sh"

echo 'cmsenv'
cmsenv

##don't setup CRAB3 unless running grid.
##Mixins.py error pops up because crab needs python 2.6, while CMSSW uses 2.7
##look up crab bootstrap script if a cmssw/crab3 compatible environment setup needed
#echo 'setting up CRAB3'
#source "$VO_CMS_SW_DIR/crab3/crab.sh"

echo 'obtaining voms proxy'
voms-proxy-init -voms cms

##alias specific to this exact cmssw setup
alias cdsrc="cd $CMSSW_BASE/src"
