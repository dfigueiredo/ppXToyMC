/****************************************************************************
 * Authors:
 *   D. Figueiredo
 ****************************************************************************/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <map>
#include <string>
#include <memory>
#include "TTree.h"
#include "Math/VectorUtil.h"

//----------------------------------------------------------------------------------------------------

class RecoParticleShow : public edm::one::EDAnalyzer<>
{
  public:
    explicit RecoParticleShow(const edm::ParameterSet&);

  private:
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void beginJob() override;
    void endJob() override;

    //RecoParticles
    edm::EDGetTokenT<std::vector<reco::PFCandidate>> RecoPartToken_;
    std::vector<const reco::PFCandidate*> recoparticlelist;

    //RecoJets
    edm::EDGetTokenT<std::vector<reco::PFJet>> RecoJetToken_;
    std::vector<const reco::PFJet*> recojetlist;

    //GenParticles
    edm::EDGetTokenT<std::vector<reco::GenParticle>> GenPartToken_;
    std::vector<const reco::GenParticle*> genparticlelist;

    //GenJets
    edm::EDGetTokenT<std::vector<reco::GenJet>> GenJetToken_;
    std::vector<const reco::GenJet*> genjetlist;


    double EBeam_;
    double MatchingMC_;
    bool DebugProtons_;
    bool DebugPF_;
    bool DebugJets_;
    double dR;

    TTree * tree_;

    std::vector<int> *recoparticle_pdgid_;
    std::vector<double> *recoparticle_energy_;
    std::vector<double> *recoparticle_pt_;
    std::vector<double> *recoparticle_eta_;
    std::vector<double> *recoparticle_phi_;
    std::vector<double> *recoparticle_px_;
    std::vector<double> *recoparticle_py_;
    std::vector<double> *recoparticle_pz_;
    std::vector<double> *recoparticle_xi_;

    std::vector<double> *recojet_energy_;
    std::vector<double> *recojet_pt_;
    std::vector<double> *recojet_eta_;
    std::vector<double> *recojet_phi_;
    std::vector<double> *recojet_px_;
    std::vector<double> *recojet_py_;
    std::vector<double> *recojet_pz_;

    std::vector<double> *genparticle_energy_;
    std::vector<double> *genparticle_pt_;
    std::vector<double> *genparticle_eta_;
    std::vector<double> *genparticle_phi_;
    std::vector<double> *genparticle_px_;
    std::vector<double> *genparticle_py_;
    std::vector<double> *genparticle_pz_;
    std::vector<double> *genparticle_xi_;
    std::vector<int> *genparticle_status_;
    std::vector<int> *genparticle_pdgid_;

    std::vector<double> *genjet_energy_;
    std::vector<double> *genjet_pt_;
    std::vector<double> *genjet_eta_;
    std::vector<double> *genjet_phi_;
    std::vector<double> *genjet_px_;
    std::vector<double> *genjet_py_;
    std::vector<double> *genjet_pz_;
    std::vector<int> *genjet_status_;
    std::vector<int> *genjet_pdgid_;

    struct orderAbsolutPz
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (fabs(vec1->pz()) > fabs(vec2->pz()));
	}
    };

    struct orderPt
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (vec1->pt() > vec2->pt());
	}
    };

};

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

RecoParticleShow::RecoParticleShow(const edm::ParameterSet& iConfig):
  RecoPartToken_        ( consumes<std::vector<reco::PFCandidate>>                                  ( iConfig.getParameter<edm::InputTag>( "RecoPartTag" ) ) ),
  RecoJetToken_         ( consumes<std::vector<reco::PFJet>>                                        ( iConfig.getParameter<edm::InputTag>( "RecoJetTag" ) ) ),
  GenPartToken_         ( consumes<std::vector<reco::GenParticle>>                                  ( iConfig.getParameter<edm::InputTag>( "GenPartTag" ) ) ),
  GenJetToken_          ( consumes<std::vector<reco::GenJet>>                                       ( iConfig.getParameter<edm::InputTag>( "GenJetTag" ) ) ),
  EBeam_   ( iConfig.getParameter<double>("EBeam")),
  MatchingMC_   ( iConfig.getParameter<bool>("MatchingMC")),
  DebugProtons_   ( iConfig.getParameter<bool>("DebugProtons")),
  DebugPF_   ( iConfig.getParameter<bool>("DebugPF")),
  DebugJets_   ( iConfig.getParameter<bool>("DebugJets"))
{
}

//----------------------------------------------------------------------------------------------------

void RecoParticleShow::analyze(const edm::Event& iEvent, const edm::EventSetup &iSetup)
{

  dR = 0;

  (*recoparticle_pdgid_).clear();
  (*recoparticle_energy_).clear();
  (*recoparticle_pt_).clear();
  (*recoparticle_eta_).clear();
  (*recoparticle_phi_).clear();
  (*recoparticle_px_).clear();
  (*recoparticle_py_).clear();
  (*recoparticle_pz_).clear();
  (*recoparticle_xi_).clear();

  recoparticlelist.clear();
  recojetlist.clear();

  (*genparticle_status_).clear();
  (*genparticle_energy_).clear();
  (*genparticle_pt_).clear();
  (*genparticle_eta_).clear();
  (*genparticle_phi_).clear();
  (*genparticle_px_).clear();
  (*genparticle_py_).clear();
  (*genparticle_pz_).clear();
  (*genparticle_xi_).clear();
  (*genparticle_pdgid_).clear();

  genparticlelist.clear();
  genjetlist.clear();

  // Get the RecoParticle collection from the event
  edm::Handle<std::vector<reco::PFCandidate>> RecoPartColl;
  iEvent.getByToken( RecoPartToken_, RecoPartColl );
  for (unsigned int i = 0; i < RecoPartColl->size(); ++i ) {
    const reco::PFCandidate* recopart = &((*RecoPartColl)[i]);
    recoparticlelist.push_back(recopart);
  }

  // Get the Jet collection from the event
  edm::Handle<std::vector<reco::PFJet> > RecoJetColl; // PAT
  iEvent.getByToken( RecoJetToken_, RecoJetColl );
  for (unsigned int i = 0; i < RecoJetColl->size(); ++i ) {
    const reco::PFJet* recojet = &((*RecoJetColl)[i]);
    recojetlist.push_back(recojet);
  }

  // Sorting PF by pz
  std::sort(recoparticlelist.begin(), recoparticlelist.end(), orderAbsolutPz());

  // Sorting Jets by pt
  std::sort(recojetlist.begin(), recojetlist.end(), orderPt());

  if(!MatchingMC_){
    for (const auto& recoparticleEvt: recoparticlelist){
      (*recoparticle_pdgid_).push_back(recoparticleEvt->pdgId());
      (*recoparticle_energy_).push_back(recoparticleEvt->energy());
      (*recoparticle_pt_).push_back(recoparticleEvt->pt());
      (*recoparticle_eta_).push_back(recoparticleEvt->eta());
      (*recoparticle_phi_).push_back(recoparticleEvt->phi());
      (*recoparticle_px_).push_back(recoparticleEvt->px());
      (*recoparticle_py_).push_back(recoparticleEvt->py());
      (*recoparticle_pz_).push_back(recoparticleEvt->pz());
      (*recoparticle_xi_).push_back(1 - abs(recoparticleEvt->pz()) / EBeam_);
      if(DebugPF_){
	std::cout << " -- Particle: \t\tpt() " << recoparticleEvt->pt() << " [GeV],  eta: "<< recoparticleEvt->eta() << ", pdgId: " << recoparticleEvt->pdgId() << std::endl;
      }
    }

    for (const auto& recojetEvt: recojetlist){
      (*recojet_energy_).push_back(recojetEvt->energy());
      (*recojet_pt_).push_back(recojetEvt->pt());
      (*recojet_eta_).push_back(recojetEvt->eta());
      (*recojet_phi_).push_back(recojetEvt->phi());
      (*recojet_px_).push_back(recojetEvt->px());
      (*recojet_py_).push_back(recojetEvt->py());
      (*recojet_pz_).push_back(recojetEvt->pz());
      if(DebugJets_){
	std::cout << " -- Jet: \t\tpt() " << recojetEvt->pt() << " [GeV], eta "<< recojetEvt->eta() << std::endl; 
      }
    }
  }

  try{

    // Get the GenParticle collection from the event
    edm::Handle<std::vector<reco::GenParticle>> GenPartColl;
    iEvent.getByToken( GenPartToken_, GenPartColl );
    for (unsigned int i = 0; i < GenPartColl->size(); ++i ) {
      const reco::GenParticle* genpart = &((*GenPartColl)[i]);
      genparticlelist.push_back(genpart);
    }

    // Get the Jet collection from the event
    edm::Handle<std::vector<reco::GenJet> > GenJetColl; // PAT
    iEvent.getByToken( GenJetToken_, GenJetColl );
    for (unsigned int i = 0; i < GenJetColl->size(); ++i ) {
      const reco::GenJet* genjet = &((*GenJetColl)[i]);
      genjetlist.push_back(genjet);
    }

    // Sorting PF by pz
    std::sort(genparticlelist.begin(), genparticlelist.end(), orderAbsolutPz());

    // Sorting Jets by pt
    std::sort(genjetlist.begin(), genjetlist.end(), orderPt());

    if(!MatchingMC_){

      for (const auto& genparticleEvt: genparticlelist){
	(*genparticle_status_).push_back(genparticleEvt->status());
	(*genparticle_pdgid_).push_back(genparticleEvt->pdgId());
	(*genparticle_energy_).push_back(genparticleEvt->energy());
	(*genparticle_pt_).push_back(genparticleEvt->pt());
	(*genparticle_eta_).push_back(genparticleEvt->eta());
	(*genparticle_phi_).push_back(genparticleEvt->phi());
	(*genparticle_px_).push_back(genparticleEvt->px());
	(*genparticle_py_).push_back(genparticleEvt->py());
	(*genparticle_pz_).push_back(genparticleEvt->pz());
	(*genparticle_xi_).push_back(1 - abs(genparticleEvt->pz()) / EBeam_);

	if(DebugProtons_){
	  if(genparticleEvt->pdgId()==2212 && genparticleEvt->status()==1 && fabs(genparticleEvt->pz())>EBeam_*0.5){
	    std::cout << " -- Proton: \t\tpt() " << genparticleEvt->pt() << " [GeV], pz " << genparticleEvt->pz() << " [GeV], eta: "<< genparticleEvt->eta() << ", status: " << genparticleEvt->status() << std::endl;
	  }
	}

	if(DebugPF_){
	  if(genparticleEvt->status()==1 && fabs(genparticleEvt->pdgId()!=2212)){
	    std::cout << " -- Particle: \t\tpt() " << genparticleEvt->pt() << " [GeV],  eta: "<< genparticleEvt->eta() << ", status: " << genparticleEvt->status() << ", pdgId: " << genparticleEvt->pdgId() << std::endl;
	  }

	  if(genparticleEvt->pdgId()==23){
	    std::cout << "Boson Z: " << std::endl;
	    std::cout << "\tStatus: " << genparticleEvt->status() << std::endl;
	    std::cout << "\tpT [GeV]: " << genparticleEvt->pt() << std::endl;
	    std::cout << "\teta [ua]: " << genparticleEvt->eta() << std::endl;
	    std::cout << "\tphi [ua]: " << genparticleEvt->phi() << std::endl;
	  }
	}
      }

      for (const auto& genjetEvt: genjetlist){
	(*genjet_status_).push_back(genjetEvt->status());
	(*genjet_pdgid_).push_back(genjetEvt->pdgId());
	(*genjet_energy_).push_back(genjetEvt->energy());
	(*genjet_pt_).push_back(genjetEvt->pt());
	(*genjet_eta_).push_back(genjetEvt->eta());
	(*genjet_phi_).push_back(genjetEvt->phi());
	(*genjet_px_).push_back(genjetEvt->px());
	(*genjet_py_).push_back(genjetEvt->py());
	(*genjet_pz_).push_back(genjetEvt->pz());
	if(DebugJets_){
	  std::cout << " -- Jet: \t\tpt() " << genjetEvt->pt() << " [GeV], eta "<< genjetEvt->eta() << ", status "<< genjetEvt->status() << std::endl; 
	}
      }

    }else{

      for(const auto& genparticleEvt: genparticlelist){
	dR = 0;
	for(const auto& recoparticleEvt: recoparticlelist){
	  //deltaR = ROOT::Math::VectorUtil::DeltaR(genparticleEvt->p4(), recoparticleEvt->p4());
	  dR = deltaR(genparticleEvt->eta(), genparticleEvt->phi(), recoparticleEvt->eta(), recoparticleEvt->phi());
	  if(dR < 0.1 ){
	    (*genparticle_status_).push_back(genparticleEvt->status());
	    (*genparticle_pdgid_).push_back(genparticleEvt->pdgId());
	    (*genparticle_energy_).push_back(genparticleEvt->energy());
	    (*genparticle_pt_).push_back(genparticleEvt->pt());
	    (*genparticle_eta_).push_back(genparticleEvt->eta());
	    (*genparticle_phi_).push_back(genparticleEvt->phi());
	    (*genparticle_px_).push_back(genparticleEvt->px());
	    (*genparticle_py_).push_back(genparticleEvt->py());
	    (*genparticle_pz_).push_back(genparticleEvt->pz());
	    (*genparticle_xi_).push_back(1 - abs(genparticleEvt->pz()) / EBeam_);
	    (*recoparticle_pdgid_).push_back(recoparticleEvt->pdgId());
	    (*recoparticle_energy_).push_back(recoparticleEvt->energy());
	    (*recoparticle_pt_).push_back(recoparticleEvt->pt());
	    (*recoparticle_eta_).push_back(recoparticleEvt->eta());
	    (*recoparticle_phi_).push_back(recoparticleEvt->phi());
	    (*recoparticle_px_).push_back(recoparticleEvt->px());
	    (*recoparticle_py_).push_back(recoparticleEvt->py());
	    (*recoparticle_pz_).push_back(recoparticleEvt->pz());
	    (*recoparticle_xi_).push_back(1 - abs(recoparticleEvt->pz()) / EBeam_);
	  } 
	}
      }

      for(const auto& genjetEvt: genjetlist){
	dR = 0;
	for(const auto& recojetEvt: recojetlist){
	  //deltaR = ROOT::Math::VectorUtil::DeltaR(genjetEvt->p4(), recojetEvt->p4());
	  dR = deltaR(genjetEvt->eta(), genjetEvt->phi(), recojetEvt->eta(), recojetEvt->phi());
	  if(dR < 0.1 ){
	    (*genjet_status_).push_back(genjetEvt->status());
	    (*genjet_energy_).push_back(genjetEvt->energy());
	    (*genjet_pt_).push_back(genjetEvt->pt());
	    (*genjet_eta_).push_back(genjetEvt->eta());
	    (*genjet_phi_).push_back(genjetEvt->phi());
	    (*genjet_px_).push_back(genjetEvt->px());
	    (*genjet_py_).push_back(genjetEvt->py());
	    (*genjet_pz_).push_back(genjetEvt->pz());
	    (*recojet_energy_).push_back(recojetEvt->energy());
	    (*recojet_pt_).push_back(recojetEvt->pt());
	    (*recojet_eta_).push_back(recojetEvt->eta());
	    (*recojet_phi_).push_back(recojetEvt->phi());
	    (*recojet_px_).push_back(recojetEvt->px());
	    (*recojet_py_).push_back(recojetEvt->py());
	    (*recojet_pz_).push_back(recojetEvt->pz());
	  } 
	}
      }
    }

  }

  catch(...){
    std::cout << "There is no GEN collection" << std::endl;
  }

  tree_->Fill();

}


//----------------------------------------------------------------------------------------------------
void RecoParticleShow::beginJob(){

  edm::Service<TFileService> fs;
  tree_=fs->make<TTree>("analyzer","analyzer");

  recoparticle_pdgid_ = new std::vector<int>;
  recoparticle_energy_ = new std::vector<double>;
  recoparticle_pt_ = new std::vector<double>;
  recoparticle_eta_ = new std::vector<double>;
  recoparticle_phi_ = new std::vector<double>;
  recoparticle_px_ = new std::vector<double>;
  recoparticle_py_ = new std::vector<double>;
  recoparticle_pz_ = new std::vector<double>;
  recoparticle_xi_ = new std::vector<double>;

  recojet_energy_ = new std::vector<double>;
  recojet_pt_ = new std::vector<double>;
  recojet_eta_ = new std::vector<double>;
  recojet_phi_ = new std::vector<double>;
  recojet_px_ = new std::vector<double>;
  recojet_py_ = new std::vector<double>;
  recojet_pz_ = new std::vector<double>;

  tree_->Branch("recoparticle_pdgid",&recoparticle_pdgid_);
  tree_->Branch("recoparticle_energy",&recoparticle_energy_);
  tree_->Branch("recoparticle_pt",&recoparticle_pt_);
  tree_->Branch("recoparticle_eta",&recoparticle_eta_);
  tree_->Branch("recoparticle_phi",&recoparticle_phi_);
  tree_->Branch("recoparticle_px",&recoparticle_px_);
  tree_->Branch("recoparticle_py",&recoparticle_py_);
  tree_->Branch("recoparticle_pz",&recoparticle_pz_);
  tree_->Branch("recoparticle_xi",&recoparticle_xi_);
  tree_->Branch("recojet_energy",&recojet_energy_);
  tree_->Branch("recojet_pt",&recojet_pt_);
  tree_->Branch("recojet_eta",&recojet_eta_);
  tree_->Branch("recojet_phi",&recojet_phi_);
  tree_->Branch("recojet_px",&recojet_px_);
  tree_->Branch("recojet_py",&recojet_py_);
  tree_->Branch("recojet_pz",&recojet_pz_);

  genparticle_status_ = new std::vector<int>;
  genparticle_pdgid_ = new std::vector<int>;
  genparticle_energy_ = new std::vector<double>;
  genparticle_pt_ = new std::vector<double>;
  genparticle_eta_ = new std::vector<double>;
  genparticle_phi_ = new std::vector<double>;
  genparticle_px_ = new std::vector<double>;
  genparticle_py_ = new std::vector<double>;
  genparticle_pz_ = new std::vector<double>;
  genparticle_xi_ = new std::vector<double>;

  genjet_status_ = new std::vector<int>;
  genjet_pdgid_ = new std::vector<int>;
  genjet_energy_ = new std::vector<double>;
  genjet_pt_ = new std::vector<double>;
  genjet_eta_ = new std::vector<double>;
  genjet_phi_ = new std::vector<double>;
  genjet_px_ = new std::vector<double>;
  genjet_py_ = new std::vector<double>;
  genjet_pz_ = new std::vector<double>;

  if(MatchingMC_){
    tree_->Branch("genparticle_status",&genparticle_status_);
    tree_->Branch("genparticle_pdgid",&genparticle_pdgid_);
    tree_->Branch("genparticle_energy",&genparticle_energy_);
    tree_->Branch("genparticle_pt",&genparticle_pt_);
    tree_->Branch("genparticle_eta",&genparticle_eta_);
    tree_->Branch("genparticle_phi",&genparticle_phi_);
    tree_->Branch("genparticle_px",&genparticle_px_);
    tree_->Branch("genparticle_py",&genparticle_py_);
    tree_->Branch("genparticle_pz",&genparticle_pz_);
    tree_->Branch("genparticle_xi",&genparticle_xi_);
    tree_->Branch("genjet_status",&genjet_status_);
    tree_->Branch("genjet_pdgid",&genjet_pdgid_);
    tree_->Branch("genjet_energy",&genjet_energy_);
    tree_->Branch("genjet_pt",&genjet_pt_);
    tree_->Branch("genjet_eta",&genjet_eta_);
    tree_->Branch("genjet_phi",&genjet_phi_);
    tree_->Branch("genjet_px",&genjet_px_);
    tree_->Branch("genjet_py",&genjet_py_);
    tree_->Branch("genjet_pz",&genjet_pz_);
  }

}


//----------------------------------------------------------------------------------------------------
void RecoParticleShow::endJob()
{

  delete recoparticle_pdgid_;
  delete recoparticle_energy_;
  delete recoparticle_pt_;
  delete recoparticle_eta_;
  delete recoparticle_phi_;
  delete recoparticle_px_;
  delete recoparticle_py_;
  delete recoparticle_pz_;
  delete recoparticle_xi_;

  delete recojet_energy_;
  delete recojet_pt_;
  delete recojet_eta_;
  delete recojet_phi_;
  delete recojet_px_;
  delete recojet_py_;
  delete recojet_pz_;

}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(RecoParticleShow);
