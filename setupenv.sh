export SCRAM_ARCH="slc6_amd64_gcc491" 

echo 'cmsenv'
cmsenv

echo 'voms initializing...'
source /osg/osg3.2/osg-client/setup.sh
voms-proxy-init -voms cms

#echo 'setting CMSSW_BASE/src alias'
alias cdsrc="cd $CMSSW_BASE/src"
