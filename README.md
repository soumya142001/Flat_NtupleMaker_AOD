# Flat_NtupleMaker_AOD
This list of files contains codes to make a flat ntuple or flat tree from an AOD file. <br>
CMSSW version : CMSSW_12_0_0 on lxplus7 (slc7_amd64_gcc*) <br>
Important files for this work are - plugins/EGamma_ntuple.cc and python/ConfFile_cfg.py 

EGamma_ntuple.cc is the analyzer code to select out branches from GsfElectron, Gsftrack, and track collection in the AOD format and put them into the flat root tree. <br>
ConfFile_cfg.py is the driver script. The AOD files on which the analyzer code should run are passed on as a parameter in the ConfFile_cfg.py code (process.source = cms.poolsource(...)) <br>

The output of the driver script is the flat root tree file (.root) <br>

The list of commands one should follow to arrive at this output .root file are: <br>
$ cd plugins <br>
$ scram b (This command will compile the .cc analyzer code and will give errors if the code has any syntax errors) <br>
$ cd ../python <br>
$ cmsRun ConfFile_cfg.py (This command will run the driver .py file. Any errors that show up most of the time are runtime errors) <br>
