// -*- C++ -*-
//
// Package:    EGamma_Tree/EGamma_ntuple
// Class:      EGamma_ntuple
//
/**\class EGamma_ntuple EGamma_ntuple.cc EGamma_Tree/EGamma_ntuple/plugins/EGamma_ntuple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Soumya Sarkar
//         Created:  Wed, 06 Sep 2023 05:56:06 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <TLorentzVector.h>

//Utilities Headers                                                                                                                                                                                      
#include <memory>
#include <string>
#include <map>
#include <iostream>
#include <cmath>
#include <fstream>
#include <TClonesArray.h>
#include <TObject.h>
#include <TChain.h>
#include <TROOT.h>
#include <vector>
#include "TLorentzVector.h"
//#include "TH1.h"                                                                                                                                                                                       
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH2.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//Gsf electrons Include Files                                                                                                                                                                            
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//Photons include files                                                                                                                                                                  
#include "DataFormats/EgammaCandidates/interface/Photon.h"
//Conversions Include file                                                                                                                                                                               
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
//DeDx hit info include file                                                                                                                                                                             
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"

//Track include files                                                                                                                                                                                    
#include "DataFormats/TrackReco/interface/Track.h"

//Calocluster and supercluster include files                                                                                                                                                            
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
//EcalRecHit class                                                                                                   
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
//GsfTrack include files                                                                                                                                                                                 
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
//Genparticle header files
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


using namespace edm;
using namespace std;


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

//Defining a Structure


using reco::TrackCollection;

class EGamma_ntuple : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit EGamma_ntuple(const edm::ParameterSet&);
  ~EGamma_ntuple();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  
  //edm::InputTag
  bool isData_;
 
  edm::EDGetTokenT<vector<reco::GsfElectron> > gsf_token_;
  edm::EDGetTokenT<vector<reco::Track> > track_token_;
  edm::EDGetTokenT<vector<reco::GsfTrack> > gsftrack_token_;
  edm::EDGetTokenT<vector<reco::GenParticle> > gen_token_;
  TTree *tree;
  
  unsigned int run_,event_,lumi_;
  
  unsigned int ngsf_;
  vector<float> gsfPx_;
  vector<float> gsfPy_;
  vector<float> gsfPz_;
  vector<float> gsfE_;
  vector<float> gsfPt_;
  vector<float> gsfEta_;
  vector<float> gsfPhi_;
  vector<int> gsfId_;
  vector<float> gsfscPixCharge_;
  vector<int> gsfseedtrackcharge_;
  vector<float> gsfseedtrackPtmode_;
  vector<float> gsfseedtrackEtamode_;
  vector<float> gsfseedtrackPhimode_;
  vector<float> gsfseedtrackPmode_;
  vector<float> gsfseedtrackPxmode_;
  vector<float> gsfseedtrackPymode_;
  vector<float> gsfseedtrackPzmode_;
  vector<float> gsfseedtrackPt_;
  vector<float> gsfseedtrackEta_;
  vector<float> gsfseedtrackPhi_;
  vector<float> gsfseedtrackP_;
  vector<float> gsfseedtrackPx_;
  vector<float> gsfseedtrackPy_;
  vector<float> gsfseedtrackPz_;
  vector<float> gsfsigmaEtaEta_;
  vector<float> gsfsigmaIetaIeta_;
  vector<float> gsfsigmaIphiIphi_;
  vector<float> gsfe1x5_;
  vector<float> gsfr9_;
  vector<float> gsfe5x5_;
  vector<float> gsfe2x5Max_;
  vector<float> gsfcorrectedEcalEnergy_;
  vector<float> gsfrawEnergy_;
  vector<float> gsfeSuperClusterOverP_;
  vector<float> gsfeSeedClusterOverP_;
  vector<float> gsfeSeedClusterOverPout_;
  vector<float> gsfeEleClusterOverPout_;
  vector<float> gsfdeltaEtaSuperClusterTrackAtVtx_;
  vector<float> gsfdeltaEtaSeedClusterTrackAtCalo_;
  vector<float> gsfdeltaEtaEleClusterTrackAtCalo_;
  vector<float> gsfdeltaPhiEleClusterTrackAtCalo_;
  vector<float> gsfdeltaPhiSuperClusterTrackAtVtx_;
  vector<float> gsfdeltaPhiSeedClusterTrackAtCalo_;
  vector<float> gsfdeltaEtaSeedClusterTrackAtVtx_;
  vector<float> gsfdxy_;
  vector<float> gsfdz_;
  vector<int> gsfseedieta_;
  vector<int> gsfseediphi_;
  vector<float> gsfseedclusterenergy_;
  vector<float> gsfseedclustereta_;
  vector<float> gsfseedclusterphi_;
  vector<Bool_t> gsfisEB_;
  vector<Bool_t> gsfisEE_;
  vector<float> gsfdr03TkSumPt_;
  vector<float> gsfdr03TkSumPtHEEP_;
  vector<float> gsfdr03EcalRecHitSumEt_;
  vector<float> gsfdr04TkSumPt_;
  vector<float> gsfdr04TkSumPtHEEP_;
  vector<float> gsfdr04EcalRecHitSumEt_;
  vector<float> gsfdr04HcalTowerSumEt_;
  vector<float> gsfdr04HcalTowerSumEtBc_;
  vector<float> gsfhcalOverEcal_;
  vector<float> gsfhcalOverEcalBc_;
  vector<float> gsffull5x5_sigmaEtaEta_;
  vector<float> gsffull5x5_sigmaIetaIeta_;
  vector<float> gsffull5x5_sigmaIphiIphi_;
  vector<float> gsffull5x5_e1x5_;
  vector<float> gsffull5x5_e2x5Max_;
  vector<float> gsffull5x5_e5x5_;
  vector<float> gsffull5x5_r9_;
  vector<float> gsffull5x5_hcalOverEcal_;
  vector<float> gsffull5x5_hcalOverEcalBc_;
  vector<Bool_t> gsfambiguous_;

  //Track Variables
  
  unsigned int ntrack_;
  vector<float> trackchi2_;
  vector<float> trackndof_;
  vector<float> tracknormalizedChi2_;
  vector<int> trackcharge_;
  vector<float> trackqoverp_;
  vector<float> tracktheta_;
  vector<float> tracklambda_;
  vector<float> trackdxy_;
  vector<float> trackdz_;
  vector<float> trackp2_;
  vector<float> trackp_;
  vector<float> trackpt2_;
  vector<float> trackpt_;
  vector<float> trackpx_;
  vector<float> trackpy_;
  vector<float> trackpz_;
  vector<float> trackphi_;
  vector<float> tracketa_;
  vector<float> trackvx_; //X-coordinate of reference point on track
  vector<float> trackvy_;
  vector<float> trackvz_;
  vector<float> trackbeta_;
  vector<float> trackqoverpError_;
  vector<float> trackptError2_;
  vector<float> trackptError_;
  vector<float> trackthetaError_;
  vector<float> tracklambdaError_;
  vector<float> tracketaError_;
  vector<float> trackphiError_;
  vector<float> trackdxyError_;
  vector<float> trackdzError_;
  vector<float> trackbetaError_;
  vector<int> tracknumberOfValidHits_;
  vector<int> tracknumberOfLostHits_;
  vector<int> trackmissingInnerHits_;
  vector<int> trackmissingOuterHits_;
  vector<float> trackvalidFraction_;
  vector<int> trackrecHitsSize_;

  //GsfTack variables

  unsigned int ngsftrack_;
  vector<float> gsftrackchargeMode_;
  vector<float> gsftrackqoverpMode_;
  vector<float> gsftrackthetaMode_;
  vector<float> gsftracklambdaMode_;
  vector<float> gsftrackpMode_;
  vector<float> gsftrackptMode_;
  vector<float> gsftrackpxMode_;
  vector<float> gsftrackpyMode_;
  vector<float> gsftrackpzMode_;
  vector<float> gsftrackphiMode_;
  vector<float> gsftracketaMode_;
  vector<float> gsftrackqoverpModeError_;
  vector<float> gsftrackptModeError_;
  vector<float> gsftrackthetaModeError_;
  vector<float> gsftracklambdaModeError_;
  vector<float> gsftracketaModeError_;
  vector<float> gsftrackphiModeError_;
  vector<float> gsftrackPt_;
  vector<float> gsftrackEta_;
  vector<float> gsftrackPhi_;
  vector<float> gsftrackP_;
  vector<float> gsftrackPx_;
  vector<float> gsftrackPy_;
  vector<float> gsftrackPz_;
  vector<float> gsftrackdxy_;
  vector<float> gsftrackdz_;
  vector<int> gsftrackcharge_;
  vector<float> gsftrackchi2_;
  vector<int> gsftrackndof_;
  vector<float> gsftracknormalizedChi2_;


  //GenParticle variables

  unsigned int ngenparticle_;
  vector<int> nmother_;
  vector<int> ndaughter_;
  vector<int> genstatus_;
  vector<int> pdgid_;
  vector<int> momid_;
  vector<vector<int>> dauid_;
  vector<float> genpT_;
  vector<float> genEta_;
  vector<float> genPhi_;
  vector<float> genM_;
  vector<float> genE_;
  vector<float> genvx_;
  vector<float> genvy_;
  vector<float> genvz_;
  

  

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EGamma_ntuple::EGamma_ntuple(const edm::ParameterSet& iConfig){

  isData_ = iConfig.getParameter<bool>("isData");

  gsf_token_ = (consumes<vector<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("ged_gsf_token")));
  track_token_ = (consumes<vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("track_token")));
  gsftrack_token_ = (consumes<vector<reco::GsfTrack> >(iConfig.getParameter<edm::InputTag>("gsftrack_token")));
  
  if(!isData_){ 
    gen_token_ = (consumes<vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("gen_token")));
  }
  
  edm::Service<TFileService> fs;
  
  //Defining Tree and tree branches
  tree = fs->make<TTree>("tree","tree");
  
  tree->Branch("run",&run_,"run/i");
  tree->Branch("event",&event_,"event/i");
  tree->Branch("lumi",&lumi_,"lumi/i");

  //Setting GsfElectron collection Branches

  tree->Branch("ngsf",&ngsf_,"ngsf/i");
  tree->Branch("gsfPx","std::vector<float>",&gsfPx_);
  tree->Branch("gsfPy","std::vector<float>",&gsfPy_);
  tree->Branch("gsfPz","std::vector<float>",&gsfPz_);
  tree->Branch("gsfE","std::vector<float>",&gsfE_);
  tree->Branch("gsfPt","std::vector<float>",&gsfPt_);
  tree->Branch("gsfEta","std::vector<float>",&gsfEta_);
  tree->Branch("gsfPhi","std::vector<float>",&gsfPhi_);
  tree->Branch("gsfId","std::vector<int>",&gsfId_);
  tree->Branch("gsfscPixCharge","std::vector<float>",&gsfscPixCharge_);
  tree->Branch("gsfseedtrackcharge","std::vector<int>",&gsfseedtrackcharge_);
  tree->Branch("gsfseedtrackPtmode","std::vector<float>",&gsfseedtrackPtmode_);
  tree->Branch("gsfseedtrackEtamode","std::vector<float>",&gsfseedtrackEtamode_);
  tree->Branch("gsfseedtrackPhimode","std::vector<float>",&gsfseedtrackPhimode_);
  tree->Branch("gsfseedtrackPmode","std::vector<float>",&gsfseedtrackPmode_);
  tree->Branch("gsfseedtrackPxmode","std::vector<float>",&gsfseedtrackPxmode_);
  tree->Branch("gsfseedtrackPymode","std::vector<float>",&gsfseedtrackPymode_);
  tree->Branch("gsfseedtrackPzmode","std::vector<float>",&gsfseedtrackPzmode_);
  tree->Branch("gsfseedtrackPt","std::vector<float>",&gsfseedtrackPt_);
  tree->Branch("gsfseedtrackEta","std::vector<float>",&gsfseedtrackEta_);
  tree->Branch("gsfseedtrackPhi","std::vector<float>",&gsfseedtrackPhi_);
  tree->Branch("gsfseedtrackP","std::vector<float>",&gsfseedtrackP_);
  tree->Branch("gsfseedtrackPx","std::vector<float>",&gsfseedtrackPx_);
  tree->Branch("gsfseedtrackPy","std::vector<float>",&gsfseedtrackPy_);
  tree->Branch("gsfseedtrackPz","std::vector<float>",&gsfseedtrackPz_);
  tree->Branch("gsfsigmaEtaEta","std::vector<float>",&gsfsigmaEtaEta_);
  tree->Branch("gsfsigmaIetaIeta","std::vector<float>",&gsfsigmaIetaIeta_);
  tree->Branch("gsfsigmaphiIphi","std::vector<float>",&gsfsigmaIphiIphi_);
  tree->Branch("gsfe1x5","std::vector<float>",&gsfe1x5_);
  tree->Branch("gsfe5x5","std::vector<float>",&gsfe5x5_);
  tree->Branch("gsfe2x5Max","std::vector<float>",&gsfe2x5Max_);
  tree->Branch("gsfr9","std::vector<float>",&gsfr9_);
  tree->Branch("gsfeSuperClusterOverP","std::vector<float>",&gsfeSuperClusterOverP_);
  tree->Branch("gsfeSeedClusterOverP","std::vector<float>",&gsfeSeedClusterOverP_);
  tree->Branch("gsfeSeedClusterOverPout","std::vector<float>",&gsfeSeedClusterOverPout_);
  tree->Branch("gsfeEleClusterOverPout","std::vector<float>",&gsfeEleClusterOverPout_);
  tree->Branch("gsfdeltaEtaSuperClusterTrackAtVtx","std::vector<float>",&gsfdeltaEtaSuperClusterTrackAtVtx_);
  tree->Branch("gsfdeltaEtaSeedClusterTrackAtCalo","std::vector<float>",&gsfdeltaEtaSeedClusterTrackAtCalo_);
  tree->Branch("gsfdeltaEtaEleClusterTrackAtCalo","std::vector<float>",&gsfdeltaEtaEleClusterTrackAtCalo_);
  tree->Branch("gsfdeltaPhiEleClusterTrackAtCalo","std::vector<float>",&gsfdeltaPhiEleClusterTrackAtCalo_);
  tree->Branch("gsfdeltaPhiSuperClusterTrackAtVtx","std::vector<float>",&gsfdeltaPhiSuperClusterTrackAtVtx_);
  tree->Branch("gsfdeltaPhiSeedClusterTrackAtCalo","std::vector<float>",&gsfdeltaPhiSeedClusterTrackAtCalo_);
  tree->Branch("gsfdeltaEtaSeedClusterTrackAtVtx","std::vector<float>",&gsfdeltaEtaSeedClusterTrackAtVtx_);
  tree->Branch("gsfrawEnergy","std::vector<float>",&gsfrawEnergy_);
  tree->Branch("gsfcorrectedEcalEnergy","std::vector<float>",&gsfcorrectedEcalEnergy_);
  tree->Branch("gsfseedieta","std::vector<int>",&gsfseedieta_);
  tree->Branch("gsfseediphi","std::vector<int>",&gsfseediphi_);
  tree->Branch("gsfseedclusterenergy","std::vector<float>",&gsfseedclusterenergy_);
  tree->Branch("gsfseedclustereta","std::vector<float>",&gsfseedclustereta_);
  tree->Branch("gsfseedclusterphi","std::vector<float>",&gsfseedclusterphi_);
  tree->Branch("gsfisEB","std::vector<Bool_t>",&gsfisEB_);
  tree->Branch("gsfisEE","std::vector<Bool_t>",&gsfisEE_);
  tree->Branch("gsfdr03TkSumPt","std::vector<float>",&gsfdr03TkSumPt_);
  tree->Branch("gsfdr03TkSumPtHEEP","std::vector<float>",&gsfdr03TkSumPtHEEP_);
  tree->Branch("gsfdr03EcalRecHitSumEt","std::vector<float>",&gsfdr03EcalRecHitSumEt_);
  tree->Branch("gsfdr04TkSumPt","std::vector<float>",&gsfdr04TkSumPt_);
  tree->Branch("gsfdr04TkSumPtHEEP","std::vector<float>",&gsfdr04TkSumPtHEEP_);
  tree->Branch("gsfdr04EcalRecHitSumEt","std::vector<float>",&gsfdr04EcalRecHitSumEt_);
  tree->Branch("gsfdr04HcalTowerSumEt","std::vector<float>",&gsfdr04HcalTowerSumEt_);
  tree->Branch("gsfdr04HcalTowerSumEtBc","std::vector<float>",&gsfdr04HcalTowerSumEtBc_);
  tree->Branch("gsfhcalOverEcal","std::vector<float>",&gsfhcalOverEcal_);
  tree->Branch("gsfhcalOverEcalBc","std::vector<float>",&gsfhcalOverEcalBc_);
  tree->Branch("gsffull5x5_sigmaEtaEta","std::vector<float>",&gsffull5x5_sigmaEtaEta_);
  tree->Branch("gsffull5x5_sigmaIetaIeta","std::vector<float>",&gsffull5x5_sigmaIetaIeta_);
  tree->Branch("gsffull5x5_sigmaIphiIphi","std::vector<float>",&gsffull5x5_sigmaIphiIphi_);
  tree->Branch("gsffull5x5_e1x5","std::vector<float>",&gsffull5x5_e1x5_);
  tree->Branch("gsffull5x5_e2x5Max","std::vector<float>",&gsffull5x5_e2x5Max_);
  tree->Branch("gsffull5x5_e5x5","std::vector<float>",&gsffull5x5_e5x5_);
  tree->Branch("gsffull5x5_r9","std::vector<float>",&gsffull5x5_r9_);
  tree->Branch("gsffull5x5_hcalOverEcal","std::vector<float>",&gsffull5x5_hcalOverEcal_);
  tree->Branch("gsffull5x5_hcalOverEcalBc","std::vector<float>",&gsffull5x5_hcalOverEcalBc_);
  tree->Branch("gsfambiguous","std::vector<Bool_t>",&gsfambiguous_);

  //Setting Track collection branches

  tree->Branch("ntrack",&ntrack_,"ntrack/i");
  tree->Branch("trackchi2","vector<float>",&trackchi2_);
  tree->Branch("trackndof","vector<float>",&trackndof_);
  tree->Branch("tracknormalizedChi2","vector<float>",&tracknormalizedChi2_);;
  tree->Branch("trackcharge","vector<int>",&trackcharge_);
  tree->Branch("trackqoverp","vector<float>",&trackqoverp_);
  tree->Branch("tracktheta","vector<float>",&tracktheta_);
  tree->Branch("tracklambda","vector<float>",&tracklambda_);
  tree->Branch("trackdxy","vector<float>",&trackdxy_);
  tree->Branch("trackdz","vector<float>",&trackdz_);
  tree->Branch("trackp2","vector<float>",&trackp2_);
  tree->Branch("trackp","vector<float>",&trackp_);
  tree->Branch("trackpt2","vector<float>",&trackpt2_);
  tree->Branch("trackpt","vector<float>",&trackpt_);
  tree->Branch("trackpx","vector<float>",&trackpx_);
  tree->Branch("trackpy","vector<float>",&trackpy_);
  tree->Branch("trackpz","vector<float>",&trackpz_);
  tree->Branch("trackphi","vector<float>",&trackphi_);
  tree->Branch("tracketa","vector<float>",&tracketa_);
  tree->Branch("trackvx","vector<float>",&trackvx_); 
  tree->Branch("trackvy","vector<float>",&trackvy_);
  tree->Branch("trackvz","vector<float>",&trackvz_);
  tree->Branch("trackbeta","vector<float>",&trackbeta_);
  tree->Branch("trackpError","vector<float>",&trackqoverpError_);
  tree->Branch("trackptError2","vector<float>",&trackptError2_);
  tree->Branch("trackptError","vector<float>",&trackptError_);
  tree->Branch("trackthetaError","vector<float>",&trackthetaError_);
  tree->Branch("tracklambdaError","vector<float>",&tracklambdaError_);
  tree->Branch("tracketaError","vector<float>",&tracketaError_);
  tree->Branch("trackphiError","vector<float>",&trackphiError_);
  tree->Branch("trackdxyError","vector<float>",&trackdxyError_);
  tree->Branch("trackdzError","vector<float>",&trackdzError_);
  tree->Branch("trackbetaError","vector<float>",&trackbetaError_);
  tree->Branch("tracknumberofValidHits","vector<int>",&tracknumberOfValidHits_);
  tree->Branch("tracknumberOfLostHits","vector<int>",&tracknumberOfLostHits_);
  tree->Branch("trackmissingInnerHits","vector<int>",&trackmissingInnerHits_);
  tree->Branch("trackmissingOuterHits","vector<int>",&trackmissingOuterHits_);
  tree->Branch("trackvalidFraction","vector<float>",&trackvalidFraction_);
  tree->Branch("trackrecHitSize","vector<int>",&trackrecHitsSize_);

  //GsfTrack variables 

  tree->Branch("ngsftrack",&ngsftrack_,"ngsftrack/i");
  tree->Branch("gsftrackchargeMode","vector<float>",&gsftrackchargeMode_);
  tree->Branch("gsftrackqoverpMode","vector<float>",&gsftrackqoverpMode_);
  tree->Branch("gsftrackthetaMode","vector<float>",&gsftrackthetaMode_);
  tree->Branch("gsftracklambdaMode","vector<float>",&gsftracklambdaMode_);
  tree->Branch("gsftrackpMode","vector<float>",&gsftrackpMode_);
  tree->Branch("gsftrackptMode","vector<float>",&gsftrackptMode_);
  tree->Branch("gsftrackpxMode","vector<float>",&gsftrackpxMode_);
  tree->Branch("gsftrackpyMode","vector<float>",&gsftrackpyMode_);
  tree->Branch("gsftrackpzMode","vector<float>",&gsftrackpzMode_);
  tree->Branch("gsftrackphiMode","vector<float>",&gsftrackphiMode_);
  tree->Branch("gsftracketaMode","vector<float>",&gsftracketaMode_);
  tree->Branch("gsftrackpModeError","vector<float>",&gsftrackqoverpModeError_);
  tree->Branch("gsftrackptModeError","vector<float>",&gsftrackptModeError_);
  tree->Branch("gsftrackthetaModeError","vector<float>",&gsftrackthetaModeError_);
  tree->Branch("gsftracklambdaModeError","vector<float>",&gsftracklambdaModeError_);
  tree->Branch("gsftracketaModeError","vector<float>",&gsftracketaModeError_);
  tree->Branch("gsftrackphiModeError","vector<float>",&gsftrackphiModeError_);
  tree->Branch("gsftrackPt"," vector<float>", &gsftrackPt_);
  tree->Branch("gsftrackEta", "vector<float>", &gsftrackEta_);
  tree->Branch("gsftrackPhi","vector<float>",&gsftrackPhi_);
  tree->Branch("gsftrackP", "vector<float>", &gsftrackP_);
  tree->Branch("gsftrackPx","vector<float>", &gsftrackPx_);
  tree->Branch("gsftrackPy","vector<float>", &gsftrackPy_);
  tree->Branch("gsftrackPz","vector<float>", &gsftrackPz_);
  tree->Branch("gsftrackdxy","vector<float>", &gsftrackdxy_);
  tree->Branch("gsftrackdz","vector<float>", &gsftrackdz_);
  tree->Branch("gsftrackcharge","vector<int>", &gsftrackcharge_);
  tree->Branch("gsftrackchi2","vector<float>", &gsftrackchi2_);
  tree->Branch("gsftrackndof","vector<int>", &gsftrackndof_);
  tree->Branch("gsftracknormalizedChi2","vector<float>", &gsftracknormalizedChi2_);
  
  //GenPartcile variable
  if(!isData_){
    tree->Branch("ngenparticle",&ngenparticle_,"ngenparticle/i");
    tree->Branch("nmother","vector<int>",&nmother_);
    tree->Branch("ndaughter","vector<int>",&ndaughter_);
    tree->Branch("genstatus","vector<int>",&genstatus_);
    tree->Branch("pdgid","vector<int>",&pdgid_);
    tree->Branch("momid","vector<int>",&momid_);
    tree->Branch("dauid","vector<vector<int>>",&dauid_);
    tree->Branch("genpT","vector<float>",&genpT_);
    tree->Branch("genEta","vector<float>",&genEta_);
    tree->Branch("genPhi","vector<float>",&genPhi_);
    tree->Branch("genM","vector<float>",&genM_);
    tree->Branch("genE","vector<float>",&genE_);
    tree->Branch("genvx","vector<float>",&genvx_);
    tree->Branch("genvy","vector<float>",&genvy_);
    tree->Branch("genvz","vector<float>",&genvz_);
  }
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

EGamma_ntuple::~EGamma_ntuple() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void EGamma_ntuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  //Event Info
  run_ = iEvent.id().run();
  event_ = iEvent.id().event();
  lumi_ = iEvent.getLuminosityBlock().luminosityBlock();

  gsfPx_.clear();
  gsfPy_.clear();
  gsfPz_.clear();
  gsfE_.clear();
  gsfPt_.clear();
  gsfEta_.clear();
  gsfPhi_.clear();
  gsfId_.clear();
  gsfscPixCharge_.clear();
  gsfsigmaEtaEta_.clear();
  gsfsigmaIetaIeta_.clear();
  gsfsigmaIphiIphi_.clear();
  gsfe1x5_.clear();
  gsfe5x5_.clear();
  gsfe2x5Max_.clear();
  gsfr9_.clear();
  gsfcorrectedEcalEnergy_.clear();
  gsfeSuperClusterOverP_.clear();
  gsfeSeedClusterOverP_.clear();
  gsfeSeedClusterOverPout_.clear();
  gsfeEleClusterOverPout_.clear();
  gsfdeltaEtaSuperClusterTrackAtVtx_.clear();
  gsfdeltaEtaSeedClusterTrackAtCalo_.clear();
  gsfdeltaEtaEleClusterTrackAtCalo_.clear();
  gsfdeltaPhiEleClusterTrackAtCalo_.clear();
  gsfdeltaPhiSuperClusterTrackAtVtx_.clear();
  gsfdeltaPhiSeedClusterTrackAtCalo_.clear();
  gsfdeltaEtaSeedClusterTrackAtVtx_.clear();
  gsfisEB_.clear();
  gsfisEE_.clear();
  gsfseedieta_.clear();
  gsfseediphi_.clear();
  gsfrawEnergy_.clear();
  gsfseedclusterenergy_.clear();
  gsfseedclustereta_.clear();
  gsfseedclusterphi_.clear();
  gsfdxy_.clear();
  gsfdz_.clear();
  gsfseedtrackcharge_.clear();
  gsfseedtrackPt_.clear();
  gsfseedtrackEta_.clear();
  gsfseedtrackPhi_.clear();
  gsfseedtrackP_.clear();
  gsfseedtrackPx_.clear();
  gsfseedtrackPz_.clear();
  gsfseedtrackPy_.clear();
  gsfseedtrackPtmode_.clear();
  gsfseedtrackEtamode_.clear();
  gsfseedtrackPhimode_.clear();
  gsfseedtrackPmode_.clear();
  gsfseedtrackPxmode_.clear();
  gsfseedtrackPzmode_.clear();
  gsfseedtrackPymode_.clear();
  gsfdr03TkSumPt_.clear();
  gsfdr03TkSumPtHEEP_.clear();
  gsfdr03EcalRecHitSumEt_.clear();
  gsfdr04TkSumPt_.clear();
  gsfdr04TkSumPtHEEP_.clear();
  gsfdr04EcalRecHitSumEt_.clear();
  gsfdr04HcalTowerSumEt_.clear();
  gsfdr04HcalTowerSumEtBc_.clear();
  gsfhcalOverEcal_.clear();
  gsfhcalOverEcalBc_.clear();
  gsffull5x5_sigmaEtaEta_.clear();
  gsffull5x5_sigmaIetaIeta_.clear();
  gsffull5x5_sigmaIphiIphi_.clear();
  gsffull5x5_e1x5_.clear();
  gsffull5x5_e2x5Max_.clear();
  gsffull5x5_e5x5_.clear();
  gsffull5x5_r9_.clear();
  gsffull5x5_hcalOverEcal_.clear();
  gsffull5x5_hcalOverEcalBc_.clear();
  gsfambiguous_.clear();

  //Track variables clearing

  trackchi2_.clear();
  trackndof_.clear();
  tracknormalizedChi2_.clear();
  trackcharge_.clear();
  trackqoverp_.clear();
  tracktheta_.clear();
  tracklambda_.clear();
  trackdxy_.clear();
  trackdz_.clear();
  trackp2_.clear();
  trackp_.clear();
  trackpt2_.clear();
  trackpt_.clear();
  trackpx_.clear();
  trackpy_.clear();
  trackpz_.clear();
  trackphi_.clear();
  tracketa_.clear();
  trackvx_.clear();
  trackvy_.clear();
  trackvz_.clear();
  trackbeta_.clear();
  trackqoverpError_.clear();
  trackptError2_.clear();
  trackptError_.clear();
  trackthetaError_.clear();
  tracklambdaError_.clear();
  tracketaError_.clear();
  trackphiError_.clear();
  trackdxyError_.clear();
  trackdzError_.clear();
  trackbetaError_.clear();
  tracknumberOfValidHits_.clear();
  tracknumberOfLostHits_.clear();
  trackmissingInnerHits_.clear();
  trackmissingOuterHits_.clear();
  trackvalidFraction_.clear();
  trackrecHitsSize_.clear();

  //Gsftrack variables clearing

  gsftrackchargeMode_.clear();
  gsftrackqoverpMode_.clear();
  gsftrackthetaMode_.clear();
  gsftracklambdaMode_.clear();
  gsftrackpMode_.clear();
  gsftrackptMode_.clear();
  gsftrackpxMode_.clear();
  gsftrackpyMode_.clear();
  gsftrackpzMode_.clear();
  gsftrackphiMode_.clear();
  gsftracketaMode_.clear();
  gsftrackqoverpModeError_.clear();
  gsftrackptModeError_.clear();
  gsftrackthetaModeError_.clear();
  gsftracklambdaModeError_.clear();
  gsftracketaModeError_.clear();
  gsftrackphiModeError_.clear();
  gsftrackPt_.clear();
  gsftrackEta_.clear();
  gsftrackPhi_.clear();
  gsftrackP_.clear();
  gsftrackPx_.clear();
  gsftrackPy_.clear();
  gsftrackPz_.clear();
  gsftrackdxy_.clear();
  gsftrackdz_.clear();
  gsftrackcharge_.clear();
  gsftrackchi2_.clear();
  gsftrackndof_.clear();
  gsftracknormalizedChi2_.clear();


  //GenParticle variables clearing

  nmother_.clear();
  ndaughter_.clear();
  genstatus_.clear();
  pdgid_.clear();
  momid_.clear();
  dauid_.clear();
  genpT_.clear();
  genEta_.clear();
  genPhi_.clear();
  genM_.clear();
  genE_.clear();
  genvx_.clear();
  genvy_.clear();
  genvz_.clear();
  
  int n = 0;

  //GsfElectron variables
  
  Handle<vector<reco::GsfElectron> > mygsf_ged;
  iEvent.getByToken(gsf_token_,mygsf_ged);
  
  for(const reco::GsfElectron &ie : *mygsf_ged){
    n++;
    //Properties directly available to GsfElectron Class
    gsfPx_.push_back(ie.px());
    gsfPy_.push_back(ie.py());
    gsfPz_.push_back(ie.pz());
    gsfE_.push_back(ie.energy());
    gsfPt_.push_back(ie.pt());
    gsfEta_.push_back(ie.eta());
    gsfPhi_.push_back(ie.phi());
    gsfId_.push_back(ie.pdgId());
    gsfscPixCharge_.push_back(ie.scPixCharge());
    gsfsigmaEtaEta_.push_back(ie.sigmaEtaEta());
    gsfsigmaIetaIeta_.push_back(ie.sigmaIetaIeta());
    gsfsigmaIphiIphi_.push_back(ie.sigmaIphiIphi());
    gsfe1x5_.push_back(ie.e1x5());
    gsfe5x5_.push_back(ie.e5x5());
    gsfe2x5Max_.push_back(ie.e2x5Max());
    gsfr9_.push_back(ie.r9());
    gsfcorrectedEcalEnergy_.push_back(ie.correctedEcalEnergy());
    gsfeSuperClusterOverP_.push_back(ie.eSuperClusterOverP());
    gsfeSeedClusterOverP_.push_back(ie.eSeedClusterOverP());
    gsfeSeedClusterOverPout_.push_back(ie.eSeedClusterOverPout());
    gsfeEleClusterOverPout_.push_back(ie.eEleClusterOverPout());
    gsfdeltaEtaSuperClusterTrackAtVtx_.push_back(ie.deltaEtaSuperClusterTrackAtVtx());
    gsfdeltaEtaSeedClusterTrackAtCalo_.push_back(ie.deltaEtaSeedClusterTrackAtCalo());
    gsfdeltaEtaEleClusterTrackAtCalo_.push_back(ie.deltaEtaEleClusterTrackAtCalo());
    gsfdeltaPhiSuperClusterTrackAtVtx_.push_back(ie.deltaPhiSuperClusterTrackAtVtx());
    gsfdeltaPhiSeedClusterTrackAtCalo_.push_back(ie.deltaPhiSeedClusterTrackAtCalo());
    gsfdeltaPhiEleClusterTrackAtCalo_.push_back(ie.deltaPhiEleClusterTrackAtCalo());
    gsfdeltaEtaSeedClusterTrackAtVtx_.push_back(ie.deltaEtaSeedClusterTrackAtVtx());
    gsfisEB_.push_back(ie.isEB());
    gsfisEE_.push_back(ie.isEE());
    gsfdr03TkSumPt_.push_back(ie.dr03TkSumPt());
    gsfdr03TkSumPtHEEP_.push_back(ie.dr03TkSumPtHEEP());
    gsfdr03EcalRecHitSumEt_.push_back(ie.dr03EcalRecHitSumEt());
    gsfdr04TkSumPt_.push_back(ie.dr04TkSumPt());
    gsfdr04TkSumPtHEEP_.push_back(ie.dr04TkSumPtHEEP());
    gsfdr04EcalRecHitSumEt_.push_back(ie.dr04EcalRecHitSumEt());
    gsfdr04HcalTowerSumEt_.push_back(ie.dr04HcalTowerSumEt());
    gsfdr04HcalTowerSumEtBc_.push_back(ie.dr04HcalTowerSumEtBc());
    gsfhcalOverEcal_.push_back(ie.hcalOverEcal());
    gsfhcalOverEcalBc_.push_back(ie.hcalOverEcalBc());
    gsffull5x5_sigmaEtaEta_.push_back(ie.full5x5_sigmaEtaEta());
    gsffull5x5_sigmaIetaIeta_.push_back(ie.full5x5_sigmaIetaIeta());
    gsffull5x5_sigmaIphiIphi_.push_back(ie.full5x5_sigmaIphiIphi());
    gsffull5x5_e1x5_.push_back(ie.full5x5_e1x5());
    gsffull5x5_e2x5Max_.push_back(ie.full5x5_e2x5Max());
    gsffull5x5_e5x5_.push_back(ie.full5x5_e5x5());
    gsffull5x5_r9_.push_back(ie.full5x5_r9());
    gsffull5x5_hcalOverEcal_.push_back(ie.full5x5_hcalOverEcal());
    gsffull5x5_hcalOverEcalBc_.push_back(ie.full5x5_hcalOverEcalBc());
    gsfambiguous_.push_back(ie.ambiguous());

    //Properties from Supercluster and GsfTrack reference to the reco::GsfElectron

    DetId id = (ie.superCluster())->seed()->seed();
    EBDetId a(id);
    gsfseedieta_.push_back(a.ieta());
    gsfseediphi_.push_back(a.iphi());
    
    gsfrawEnergy_.push_back((ie.superCluster())->rawEnergy());
    gsfseedclusterenergy_.push_back((ie.superCluster())->seed()->energy());
    gsfseedclustereta_.push_back((ie.superCluster())->seed()->eta());
    gsfseedclusterphi_.push_back((ie.superCluster())->seed()->phi());

    gsfdxy_.push_back((ie.gsfTrack())->dxy());
    gsfdz_.push_back((ie.gsfTrack())->dz());
    gsfseedtrackcharge_.push_back((ie.gsfTrack())->charge());
    gsfseedtrackPtmode_.push_back((ie.gsfTrack())->ptMode());
    gsfseedtrackEtamode_.push_back((ie.gsfTrack())->etaMode());
    gsfseedtrackPhimode_.push_back((ie.gsfTrack())->phiMode());
    gsfseedtrackPmode_.push_back((ie.gsfTrack())->pMode());
    gsfseedtrackPxmode_.push_back((ie.gsfTrack())->pxMode());
    gsfseedtrackPymode_.push_back((ie.gsfTrack())->pyMode());
    gsfseedtrackPzmode_.push_back((ie.gsfTrack())->pzMode());
    gsfseedtrackPt_.push_back((ie.gsfTrack())->pt());
    gsfseedtrackEta_.push_back((ie.gsfTrack())->eta());
    gsfseedtrackPhi_.push_back((ie.gsfTrack())->phi());
    gsfseedtrackP_.push_back((ie.gsfTrack())->p());
    gsfseedtrackPx_.push_back((ie.gsfTrack())->px());
    gsfseedtrackPy_.push_back((ie.gsfTrack())->py());
    gsfseedtrackPz_.push_back((ie.gsfTrack())->pz());

  }

  ngsf_ = n;

  //Track Variables
  Handle<vector<reco::Track> > mytrack;
  iEvent.getByToken(track_token_,mytrack);
  int  nt=0;
  for(const reco::Track &ie : *mytrack){
    nt++;
    trackchi2_.push_back(ie.chi2());
    trackndof_.push_back(ie.ndof());
    tracknormalizedChi2_.push_back(ie.normalizedChi2());
    trackcharge_.push_back(ie.charge());
    trackqoverp_.push_back(ie.qoverp());
    tracktheta_.push_back(ie.theta());
    tracklambda_.push_back(ie.lambda());
    trackdxy_.push_back(ie.dxy());
    trackdz_.push_back(ie.dz());
    trackp2_.push_back(ie.p2());
    trackp_.push_back(ie.p());
    trackpt2_.push_back(ie.pt2());
    trackpt_.push_back(ie.pt());
    trackpx_.push_back(ie.px());
    trackpy_.push_back(ie.py());
    trackpz_.push_back(ie.pz());
    trackphi_.push_back(ie.phi());
    tracketa_.push_back(ie.eta());
    trackvx_.push_back(ie.vx());
    trackvy_.push_back(ie.vy());
    trackvz_.push_back(ie.vz());
    trackbeta_.push_back(ie.beta());
    trackqoverpError_.push_back(ie.qoverpError());
    trackptError2_.push_back(ie.ptError2());
    trackptError_.push_back(ie.ptError());
    trackthetaError_.push_back(ie.thetaError());
    tracklambdaError_.push_back(ie.lambdaError());
    tracketaError_.push_back(ie.etaError());
    trackphiError_.push_back(ie.phiError());
    trackdxyError_.push_back(ie.dxyError());
    trackdzError_.push_back(ie.dzError());
    trackbetaError_.push_back(ie.betaError());
    tracknumberOfValidHits_.push_back(ie.numberOfValidHits());
    tracknumberOfLostHits_.push_back(ie.numberOfLostHits());
    trackmissingInnerHits_.push_back(ie.missingInnerHits());
    trackmissingOuterHits_.push_back(ie.missingOuterHits());
    trackvalidFraction_.push_back(ie.validFraction());
    //trackrecHitsSize_.push_back(ie.recHitsSize());
  
  }

  ntrack_ = nt;

  //GsfTrack Variables
  Handle<vector<reco::GsfTrack> > mygsftrack;
  iEvent.getByToken(gsftrack_token_,mygsftrack);

  int ngt=0;
  
  for(const reco::GsfTrack &ie : *mygsftrack){
    ngt++;
    
    gsftrackchargeMode_.push_back(ie.chargeMode());
    gsftrackqoverpMode_.push_back(ie.qoverpMode());
    gsftrackthetaMode_.push_back(ie.thetaMode());
    gsftracklambdaMode_.push_back(ie.lambdaMode());
    gsftrackpMode_.push_back(ie.pMode());
    gsftrackptMode_.push_back(ie.ptMode());
    gsftrackpxMode_.push_back(ie.pxMode());
    gsftrackpyMode_.push_back(ie.pyMode());
    gsftrackpzMode_.push_back(ie.pzMode());
    gsftrackphiMode_.push_back(ie.phiMode());
    gsftracketaMode_.push_back(ie.etaMode());
    gsftrackqoverpModeError_.push_back(ie.qoverpModeError());
    gsftrackptModeError_.push_back(ie.ptModeError());
    gsftrackthetaModeError_.push_back(ie.thetaModeError());
    gsftracklambdaModeError_.push_back(ie.lambdaModeError());
    gsftracketaModeError_.push_back(ie.etaModeError());
    gsftrackphiModeError_.push_back(ie.phiModeError());
    gsftrackPt_.push_back(ie.pt());
    gsftrackEta_.push_back(ie.eta());
    gsftrackPhi_.push_back(ie.phi());
    gsftrackP_.push_back(ie.p());
    gsftrackPx_.push_back(ie.px());
    gsftrackPy_.push_back(ie.py());
    gsftrackPz_.push_back(ie.pz());
    gsftrackdxy_.push_back(ie.dxy());
    gsftrackdz_.push_back(ie.dz());
    gsftrackcharge_.push_back(ie.charge());
    gsftrackchi2_.push_back(ie.chi2());
    gsftrackndof_.push_back(ie.ndof());
    gsftracknormalizedChi2_.push_back(ie.normalizedChi2());
  }
  
  ngsftrack_ = ngt;

  //GenParticle Variables
  if(!isData_){
  Handle<vector<reco::GenParticle> > mygenparticle;
  iEvent.getByToken(gen_token_,mygenparticle);
  
  int ngen=0;
  
  for(const reco::GenParticle &ie : *mygenparticle){
    ngen++;
   
    nmother_.push_back(ie.numberOfMothers());
    ndaughter_.push_back(ie.numberOfDaughters());
    genstatus_.push_back(ie.status());
    pdgid_.push_back(ie.pdgId());
    genpT_.push_back(ie.pt());
    genEta_.push_back(ie.eta());
    genPhi_.push_back(ie.phi());
    genM_.push_back(ie.mass());
    genE_.push_back(ie.energy());
    genvx_.push_back(ie.vx());
    genvy_.push_back(ie.vy());
    genvz_.push_back(ie.vz());

    //Filling Genparticle mother information (first different mother)
    int id = ie.pdgId();
    int nmom = ie.numberOfMothers();
    int mid = 0;
    if(nmom>0){
      const reco::Candidate *mom = ie.mother(0);
      mid = mom->pdgId();
      while(mid == id){
	int ngrandmom = mom->numberOfMothers();
	if(ngrandmom>0){
	  const reco::Candidate *grandmom = mom->mother(0);
	  int grandmomid = grandmom->pdgId();
	  mom = grandmom;
	  id = mid;
	  mid = grandmomid;
	}
      }
    }
    momid_.push_back(mid);
    
    vector<int> dauid;
    dauid.clear();

    for(size_t i=0;i<ie.numberOfDaughters();i++) dauid.push_back(ie.daughter(i)->pdgId());
   
    dauid_.push_back(dauid);
  }

  ngenparticle_ = ngen;
  
  }
  
  //Do only once
  tree->Fill();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void EGamma_ntuple::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void EGamma_ntuple::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EGamma_ntuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EGamma_ntuple);
