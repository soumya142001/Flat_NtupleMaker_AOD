import FWCore.ParameterSet.Config as cms
import sys

#linenumber = sys.argv[2]
#linenumber = int(linenumber)

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

isData = False

with open("/eos/user/s/sosarkar/condor_Jobs/Ztoee_Samples/AOD_output_v0/Ztoee_pT40to300_AODfilename.txt") as file:
    list0=file.readlines()

with open("/eos/user/s/sosarkar/condor_Jobs/AOD_output_v9/JPsiToee_pT40to250_AODfilename.txt") as file:
    list1=file.readlines()
with open("/eos/user/s/sosarkar/condor_Jobs/AOD_output_v10/JPsiToee_pT250to500_AODfilename.txt") as file:
    list2=file.readlines()
with open("/eos/user/s/sosarkar/condor_Jobs/QCD_samples/AOD_output_v0/QCD_AOD_samples_v0.txt") as file:
    list3=file.readlines()
with open("/eos/user/s/sosarkar/condor_Jobs/QCD_samples/AOD_output_v1/QCD_AOD_samples_v1.txt") as file:
    list4=file.readlines()
with open("/eos/user/s/sosarkar/condor_Jobs/QCD_samples/AOD_output_v2/QCD_AOD_samples_v2.txt") as file:
    list5=file.readlines()
with open("/eos/user/s/sosarkar/condor_Jobs/QCD_samples/AOD_output_v3/QCD_AOD_samples_v3.txt") as file:
    list6=file.readlines()
with open("/eos/user/s/sosarkar/condor_Jobs/JPsitoee_pT40to250/AOD_output_v0/JPsitoee_pT40to250_filename.txt") as file:
    list7=file.readlines()
with open("/eos/user/s/sosarkar/condor_Jobs/JPsitoee_pT250to500/AOD_output_v0/JPsitoee_pT250to500_filename.txt") as file:
    list8=file.readlines()
with open("/eos/user/s/sosarkar/condor_Jobs/Ztoee_pT40to300/AOD_output_v0/Ztoee_pT40to300_filename.txt") as file:
    list9=file.readlines()
with open("/eos/user/s/sosarkar/condor_Jobs/Ztoee_pT40to300/AOD_output_v1/Ztoee_pT40to300_filename.txt") as file:
    list10=file.readlines()
with open("/eos/user/s/sosarkar/condor_Jobs/JPsitoee_pT40to500/AOD_output_v1/JPsitoee_pT40to500_filelist.txt") as file:
    list11=file.readlines()
with open("/eos/user/s/sosarkar/condor_Jobs/JPsitoee_pT40to500/AOD_output_v0/JPsitoee_pT40to500_filename.txt") as file:
    list12=file.readlines() 
with open("/eos/user/s/sosarkar/condor_Jobs/Ztoee_pT40to300/AOD_output_v2/Ztoee_pT40to300_filename.txt") as file:
    list13=file.readlines() 
with open("/eos/user/s/sosarkar/HNL_samples/HNL_filename.txt") as file:
    list14=file.readlines()


filelist0 = []
for i in list0:
    filelist0.append('file:'+i)

filelist1 = []
for i in list1:
    filelist1.append('file:'+i)

filelist2 = []
for i in list2:
    filelist2.append('file:'+i)

filelist3 = []
for i  in list3:
    filelist3.append('file:'+i)
filelist4 = []
for i in list4:
    filelist4.append('file:'+i)

filelist5 = []
for i in list5:
    filelist5.append('file:'+i)

filelist6 = []
for i in list6:
    filelist6.append('file:'+i)

filelist7 = []
for i  in list7:
    filelist7.append('file:'+i)

filelist8 = []
for i in list8:
    filelist8.append('file:'+i)

filelist9 = []
for i in list9:
    filelist9.append('file:'+i)

filelist10 = []
for i in list10:
    filelist10.append('file:'+i)

filelist11 = []
for i in list11:
    filelist11.append('file:'+i)

filelist12 = []
for i in list12:
    filelist12.append('file:'+i)

filelist13 = []
for i in list13:
    filelist13.append('file:'+i)

filelist14 = []
for i in list14:
    filelist14.append('file:'+i)

Ztoee_pT40to300_files = filelist13+filelist10+filelist0+filelist9
JPsitoee_pT40to500_files = filelist12+filelist2+list1+filelist7+filelist8 #+filelist11 (not added ad no gen particle branch is there)
QCD_pT150to200_files = filelist6+filelist5+filelist4+filelist3
EGamma_datafile = 'file:/eos/home-s/sosarkar/EGamma_Run2018A_UL2018_v1_data_AOD.root'
HNL_M2_datafile = filelist14
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
                                #  HNL_M2_datafile
                                # EGamma_datafile
                                # Ztoee_pT40to300_files
                                # JPsitoee_pT40to500_files
                                 QCD_pT150to200_files
                                )
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('QCDpT150to200_ntuple.root')                                                                                 
                                   )

process.demo = cms.EDAnalyzer('EGamma_ntuple',
   isData = cms.bool(isData),
   ged_gsf_token = cms.InputTag('gedGsfElectrons'),
   track_token = cms.InputTag('generalTracks'),
   gsftrack_token = cms.InputTag('electronGsfTracks'),
   gen_token = cms.InputTag('genParticles'),
   SkipEvent = cms.untracked.vstring('ProductNotFound')
                              )

process.p = cms.Path(process.demo)
