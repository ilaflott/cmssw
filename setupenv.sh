
echo 'sourcing cms software'
export VO_CMS_SW_DIR="/cvmfs/cms.cern.ch" #for access to CERN CMSSW repositories
source "$VO_CMS_SW_DIR/cmsset_default.sh"
#source /osg/app/cmssoft/cms/cmsset_default.sh

echo 'cmsenv'
cmsenv

echo 'voms initializing...'
source /osg/osg3.2/osg-client/setup.sh
voms-proxy-init -voms cms