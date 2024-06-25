import CRABClient
from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
config = Configuration()

config.section_("General")
config.General.requestName      = 'EGamma_Run2_2018A_UL2018_Jan3_v0'
config.General.transferLogs     = True

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName       = 'Analysis'
config.JobType.psetName         = '/afs/cern.ch/user/s/sosarkar/MYDEMOANALYZER/CMSSW_12_0_0/src/EGamma_Tree/EGamma_ntuple/python/ConfFile_cfg.py'
#config.JobType.sendPythonFolder = True                                                                                                                                                                 
#config.JobType.numCores = 1                                                                                                                                                                            

config.section_("Data")
config.Data.inputDataset        ='/EGamma/Run2018A-15Feb2022_UL2018-v1/AOD'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
config.Data.runRange = '315257-316995'
config.Data.inputDBS            = 'global'
#config.Data.allowNonValidInputDataset = True                                                                                                                                                            
config.Data.splitting           = 'FileBased'
config.Data.totalUnits          = 1000                                                                                                                                                                  
config.Data.unitsPerJob         = 1
config.Data.outLFNDirBase       = '/store/user/sosarkar/EGamma_Run2_2018A_UL2018_Ntuples'
config.Data.publication = False
config.Data.outputDatasetTag    = 'EGamma_Run2_2018A_UL2018_Jan03_NTuples_v0'

config.section_("Site")
config.Site.storageSite         = "T3_CH_CERNBOX"
#config.Site.blacklist           = ['T2_IN_TIFR']                                                                                                                                                        
#config.Site.whitelist           = ['T2_TR_METU']                                                                                                                                                        

