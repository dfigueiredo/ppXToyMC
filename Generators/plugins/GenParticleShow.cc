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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include <map>
#include <string>
#include <memory>
#include "TTree.h"


//----------------------------------------------------------------------------------------------------

class GenParticleShow : public edm::one::EDAnalyzer<>
{
  public:
    explicit GenParticleShow(const edm::ParameterSet&);

  private:
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void beginJob() override;
    void endJob() override;

    //GenParticles
    edm::EDGetTokenT<std::vector<reco::GenParticle>> GenPartToken_;
    std::vector<const reco::GenParticle*> genparticlelist;

    //GenJets
    edm::EDGetTokenT<edm::View<reco::GenJet>> GenJetToken_;
    std::vector<const reco::GenJet*> genjetlist;

    double PzProtonThreshold_;
    double EBeam_;
    bool DebugProtons_;
    bool DebugPF_;
    bool DebugJets_;

    TTree * tree_;

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

GenParticleShow::GenParticleShow(const edm::ParameterSet& iConfig):
  GenPartToken_        ( consumes<std::vector<reco::GenParticle>>                    ( iConfig.getParameter<edm::InputTag>( "GenPartTag" ) ) ),
  GenJetToken_        ( consumes<edm::View<reco::GenJet>>                    ( iConfig.getParameter<edm::InputTag>( "GenJetTag" ) ) ),
  EBeam_   ( iConfig.getParameter<double>("EBeam")),
  DebugProtons_   ( iConfig.getParameter<bool>("DebugProtons")),
  DebugPF_   ( iConfig.getParameter<bool>("DebugPF")),
  DebugJets_   ( iConfig.getParameter<bool>("DebugJets"))
{
}

//----------------------------------------------------------------------------------------------------

void GenParticleShow::analyze(const edm::Event& iEvent, const edm::EventSetup &iSetup)
{

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

  // Get the GenParticle collection from the event
  edm::Handle<std::vector<reco::GenParticle>> GenPartColl;
  iEvent.getByToken( GenPartToken_, GenPartColl );
  for (unsigned int i = 0; i < GenPartColl->size(); ++i ) {
    const reco::GenParticle* genpart = &((*GenPartColl)[i]);
    genparticlelist.push_back(genpart);
  }

  // Get the Jet collection from the event
  edm::Handle<edm::View<reco::GenJet> > GenJetColl; // PAT
  iEvent.getByToken( GenJetToken_, GenJetColl );
  for (unsigned int i = 0; i < GenJetColl->size(); ++i ) {
    const reco::GenJet* genjet = &((*GenJetColl)[i]);
    genjetlist.push_back(genjet);
  }

  // Sorting PF by pz
  std::sort(genparticlelist.begin(), genparticlelist.end(), orderAbsolutPz());

  // Sorting Jets by pt
  std::sort(genjetlist.begin(), genjetlist.end(), orderPt());

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
	//if(fabs(genparticleEvt->pdgId()!=2212)){
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

    tree_->Fill();

  }


  //----------------------------------------------------------------------------------------------------
  void GenParticleShow::beginJob(){

    edm::Service<TFileService> fs;
    tree_=fs->make<TTree>("analyzer","analyzer");

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


  //----------------------------------------------------------------------------------------------------
  void GenParticleShow::endJob()
  {

    delete genparticle_status_;
    delete genparticle_pdgid_;
    delete genparticle_energy_;
    delete genparticle_pt_;
    delete genparticle_eta_;
    delete genparticle_phi_;
    delete genparticle_px_;
    delete genparticle_py_;
    delete genparticle_pz_;
    delete genparticle_xi_;

    delete genjet_status_;
    delete genjet_pdgid_;
    delete genjet_energy_;
    delete genjet_pt_;
    delete genjet_eta_;
    delete genjet_phi_;
    delete genjet_px_;
    delete genjet_py_;
    delete genjet_pz_;

  }

  //----------------------------------------------------------------------------------------------------

  DEFINE_FWK_MODULE(GenParticleShow);
