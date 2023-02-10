// Authors: PPS X team
// Toy MC for pp-> pSp : S->X+H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGamma.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandChiSquare.h"
#include "CLHEP/Random/RandBreitWigner.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LesHouches.h"
#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEXMLStringProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHEReader.h"
#include "ppXToyMC/Generators/plugins/particle_id.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include <cmath>


//----------------------------------------------------------------------------------------------------

class PPSToyMcLHE_XH : public edm::one::EDProducer<edm::BeginRunProducer, edm::EndRunProducer>
{

  public:
    PPSToyMcLHE_XH(const edm::ParameterSet &);
    virtual ~PPSToyMcLHE_XH();

  private:
    virtual void produce(edm::Event & e, const edm::EventSetup& es) override;
    void beginRunProduce(edm::Run &run, edm::EventSetup const &es) override;
    void endRunProduce(edm::Run &, edm::EventSetup const &) override;
    void fillEventLHE(lhef::HEPEUP &outlhe, CLHEP::HepLorentzVector &particle, int pdgid, double status, double spin);

    LHERunInfoProduct::Header SLHA();

    // input parameters
    unsigned int verbosity;
    bool decayX;
    bool decayHToBottoms;
    const double m_S_mean;  // mass of the S scalar particle, GeV
    const double m_S_gamma; // width of the S particle, GeV
    const double m_H_mean;   // mass of the H particle, GeV
    const double m_H_gamma;  // width of the H particle, GeV
    const double m_X_mean;   // mass of the X particle, GeV
    const double m_X_gamma;  // width of the X particle, GeV
    const double m_f; 	    //mass of the fermions, GeV
    const double m_b;	    // mass of the bottom quark, GeV
    const double p_beam;    // beam momentum, GeV
    const double p_z_LAB_2p_min; // minimum of p_z of the 2-proton system in the LAB frame, GeV
    const double p_z_LAB_2p_max; // maximum of p_z of the 2-proton system in the LAB frame, GeV
    const double p_T_H_min; // minimum value of H's pT in the LAB frame, GeV
    std::string lheOutputFile; // LHE filename. If empty, this tool will not produce the LHE XML file.

    std::ofstream file;
    bool writeLHE;

};

//----------------------------------------------------------------------------------------------------

using namespace edm;
using namespace std;

//----------------------------------------------------------------------------------------------------

PPSToyMcLHE_XH::PPSToyMcLHE_XH(const edm::ParameterSet& pset) :

  verbosity(pset.getUntrackedParameter<unsigned int>("verbosity", 0)),
  decayX(pset.getParameter<bool>("decayX")),
  decayHToBottoms(pset.getParameter<bool>("decayHToBottoms")),
  m_S_mean(pset.getParameter<double>("m_S_mean")),
  m_S_gamma(pset.getParameter<double>("m_S_gamma")),
  m_H_mean(pset.getParameter<double>("m_H_mean")),
  m_H_gamma(pset.getParameter<double>("m_H_gamma")),
  m_X_mean(pset.getParameter<double>("m_X_mean")),
  m_X_gamma(pset.getParameter<double>("m_X_gamma")),
  m_f(pset.getParameter<double>("m_f")),
  m_b(pset.getParameter<double>("m_b")),
  p_beam(pset.getParameter<double>("p_beam")),
  p_z_LAB_2p_min(pset.getParameter<double>("p_z_LAB_2p_min")),
  p_z_LAB_2p_max(pset.getParameter<double>("p_z_LAB_2p_max")),
  p_T_H_min(pset.getParameter<double>("p_T_H_min")),
  lheOutputFile(pset.getUntrackedParameter<std::string>("lheOutputFile"))
{

  produces<HepMCProduct>("unsmeared");
  produces<LHEEventProduct>("source");
  produces<LHERunInfoProduct, edm::Transition::BeginRun>();

  writeLHE = false;
  if (!lheOutputFile.empty()) {
    writeLHE = true;
    file.open(lheOutputFile, std::fstream::out | std::fstream::trunc);
  }

}


void PPSToyMcLHE_XH::beginRunProduce(edm::Run &run, edm::EventSetup const &) {

  // Fill HEPRUP common block and store in edm::Run
  lhef::HEPRUP heprup;

  // Set the number of processes per events
  heprup.resize(1);

  // Beam Particle ID
  heprup.IDBMUP.first = 2212;
  heprup.IDBMUP.second = 2212;

  // Beam Particle Energy
  heprup.EBMUP.first = p_beam;
  heprup.EBMUP.second = p_beam;

  // Beam Particle PDF
  heprup.PDFGUP.first = -1;
  heprup.PDFGUP.second = -1;

  // PDF set ID for the process. Useless for PYTHIA.
  heprup.PDFSUP.first = -1;
  heprup.PDFSUP.second = -1;

  heprup.IDWTUP = 3;

  // Metadata for the first production (excluding its decays. I.e pp->pSp):
  heprup.XSECUP[0] = 1.; 
  heprup.XERRUP[0] = 0;
  heprup.XMAXUP[0] = 1;
  heprup.LPRUP[0] = 1;

  std::unique_ptr<LHERunInfoProduct> runInfo(new LHERunInfoProduct(heprup));
  runInfo->addHeader(SLHA());

  if (writeLHE) std::copy(runInfo->begin(), runInfo->end(), std::ostream_iterator<std::string>(file));
  run.put(std::move(runInfo));

}

void PPSToyMcLHE_XH::endRunProduce(edm::Run &run, edm::EventSetup const &es) {
  if(writeLHE){
    file << LHERunInfoProduct::endOfFile();
    file.close();
  }
}



//----------------------------------------------------------------------------------------------------

void PPSToyMcLHE_XH::fillEventLHE(lhef::HEPEUP &outlhe, 
    CLHEP::HepLorentzVector &particle, 
    int pdgid, 
    double status, double spin) {

  int index = outlhe.NUP;
  outlhe.resize(outlhe.NUP + 1);

  outlhe.IDPRUP = 1;
  outlhe.XWGTUP = 2.8167496e-08;
  outlhe.SCALUP = 1; 
  outlhe.AQEDUP = 7.29929980e-03;
  outlhe.AQCDUP = 9.87552740e-02;

  outlhe.PUP[index][0] = particle.px();
  outlhe.PUP[index][1] = particle.py();
  outlhe.PUP[index][2] = particle.pz();
  outlhe.PUP[index][3] = particle.e();

  if(pdgid==2212 || pdgid==-2212){
    outlhe.PUP[index][4] = 9.38270000000e-01;
    outlhe.SPINUP[index] = spin; // added later
  }else{
    outlhe.PUP[index][4] = particle.m();
  }
  outlhe.IDUP[index] = pdgid;

  if (pdgid == 5) {
    outlhe.ICOLUP[index].first = 501;
    outlhe.ICOLUP[index].second = 0;
    outlhe.MOTHUP[index].first = 9;
    outlhe.MOTHUP[index].second = 9;
  }
  else if (pdgid == -5) {
    outlhe.ICOLUP[index].first = 0;
    outlhe.ICOLUP[index].second = 501; 
    outlhe.MOTHUP[index].first = 9;
    outlhe.MOTHUP[index].second = 9;
  }
  else if (pdgid == 51) { // DM SPIN 0 Particle
    outlhe.ICOLUP[index].first = 0;
    outlhe.ICOLUP[index].second = 0;
    outlhe.MOTHUP[index].first = 1;
    outlhe.MOTHUP[index].second = 2;
  }
  else if (std::abs(pdgid) == 25 || std::abs(pdgid) == particleId_X) {
    outlhe.ICOLUP[index].first = 0;
    outlhe.ICOLUP[index].second = 0;
    outlhe.MOTHUP[index].first = 5;
    outlhe.MOTHUP[index].second = 5;
  }else if(std::abs(pdgid) == 14){
    outlhe.ICOLUP[index].first = 0;
    outlhe.ICOLUP[index].second = 0;
    outlhe.MOTHUP[index].first = 6;
    outlhe.MOTHUP[index].second = 6;
  }else{
    outlhe.ICOLUP[index].first = 0;
    outlhe.ICOLUP[index].second = 0;
    outlhe.MOTHUP[index].first = 0;
    outlhe.MOTHUP[index].second = 0;
  }
  outlhe.ISTUP[index] = status;
  return;
}

void PPSToyMcLHE_XH::produce(edm::Event &e, const edm::EventSetup& es)
{

  lhef::HEPEUP hepeup;

  CLHEP::HepLorentzVector zero(0, 0, 0, 0);

  if (verbosity)
    printf("\n>> PPSToyMcLHE_XH::produce > event %llu\n", e.id().event());

  edm::Service<edm::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine* engine = &rng->getEngine(e.streamID());

  HepMC::GenEvent *fEvt = new HepMC::GenEvent();
  fEvt->set_event_number(e.id().event());

  HepMC::GenVertex *vtx = new HepMC::GenVertex(HepMC::FourVector(0., 0., 0., 0.));
  fEvt->add_vertex(vtx);

  double m_H = -1., m_S= -1., m_X = -1.;
  double px_1 = -1., py_1 = -1.;
  double px_2 = -1., py_2 = -1.;

  // 4-momenta of the outgoing particles in the LAB frame
  CLHEP::HepLorentzVector momentum_H;
  CLHEP::HepLorentzVector momentum_X;
  CLHEP::HepLorentzVector momentum_S;

  CLHEP::HepLorentzVector momentum_p1, momentum_p1i;
  CLHEP::HepLorentzVector momentum_p2, momentum_p2i;

  // Try to generate event fullfilling all criteria
  bool generationOK = false;
  for (unsigned int n_attempt = 0; n_attempt < 10000; ++n_attempt)
  {
    m_H = CLHEP::RandBreitWigner::shoot(engine, m_H_mean, m_H_gamma);
    m_S = CLHEP::RandBreitWigner::shoot(engine, m_S_mean, m_S_gamma);
    m_X = CLHEP::RandBreitWigner::shoot(engine, m_X_mean, m_X_gamma);

    px_1 = CLHEP::RandGauss::shoot(engine, 0., 0.1);
    py_1 = CLHEP::RandGauss::shoot(engine, 0., 0.1);
    px_2 = CLHEP::RandGauss::shoot(engine, 0., 0.1);
    py_2 = CLHEP::RandGauss::shoot(engine, 0., 0.1);

    if (m_H < 0. || m_S < m_H + m_X) continue;

    // Generate p_z of the 2-proton system in the LAB frame
    const double p_z_LAB_2p = CLHEP::RandFlat::shoot(engine, p_z_LAB_2p_min, p_z_LAB_2p_max);

    // Generate spherical angles in the CMS frame of the X-Z system
    const double cos_theta_c = 2. * CLHEP::RandFlat::shoot(engine) - 1.;
    const double sin_theta_c = sqrt(1. - cos_theta_c * cos_theta_c);
    const double phi_c = CLHEP::RandFlat::shoot(engine) * 2. * M_PI;

    // Determine xi's of the protons
    // proton 1: positive z momentum component
    const double xi2 = (p_z_LAB_2p + sqrt(p_z_LAB_2p*p_z_LAB_2p + m_S*m_S)) / (2. * p_beam);
    const double xi1 = m_S * m_S / (4. * p_beam * p_beam * xi2);

    if (verbosity)
    {
      printf("  m_H = %.1f\n", m_H);
      printf("  m_S = %.1f\n", m_S);
      printf("  p_z_LAB_2p = %.1f\n", p_z_LAB_2p);
      printf("  xi1 = %.3f, xi2 = %.3f\n", xi1, xi2);
      printf("  p_beam * (xi2 - xi1) = %.1f\n", p_beam * (xi2 - xi1));
    }

    // Determine momenta of the X and S particles in the CMS frame of the X-S system
    const double p_c_sq = pow(m_S*m_S - m_X*m_X - m_H*m_H, 2.) / (4. * m_S * m_S) - m_X*m_X * m_H*m_H / (m_S*m_S);
    const double p_c = sqrt(p_c_sq);

    if (verbosity)
      printf("  p_c = %.3f\n", p_c);

    CLHEP::HepLorentzVector momentum_X_CMS(
	+ p_c * sin_theta_c * cos(phi_c),
	+ p_c * sin_theta_c * sin(phi_c),
	+ p_c * cos_theta_c,
	sqrt(p_c*p_c + m_X*m_X)
	);

    CLHEP::HepLorentzVector momentum_H_CMS(
	- p_c * sin_theta_c * cos(phi_c),
	- p_c * sin_theta_c * sin(phi_c),
	- p_c * cos_theta_c,
	sqrt(p_c*p_c + m_H*m_H)
	);

    CLHEP::HepLorentzVector momentum_S_CMS(
	- p_c * sin_theta_c * cos(phi_c),
	- p_c * sin_theta_c * sin(phi_c),
	- p_c * cos_theta_c,
	sqrt(p_c*p_c + m_S*m_S)
	);

    // Determine boost from X-H CMS frame to the LAB frame
    const double beta = (xi1 - xi2) / (xi1 + xi2);
    const CLHEP::Hep3Vector betaVector(0., 0., beta);

    if (verbosity)
      printf("  beta = %.3f\n", beta);

    // Determine 4-momenta of the outgoing particles in the LAB frame
    momentum_H = CLHEP::boostOf(momentum_H_CMS, betaVector);
    momentum_X = CLHEP::boostOf(momentum_X_CMS, betaVector);
    momentum_S = CLHEP::boostOf(momentum_S_CMS, betaVector);

    momentum_p1i = CLHEP::HepLorentzVector(0., 0., +p_beam, p_beam);
    momentum_p2i = CLHEP::HepLorentzVector(0., 0., -p_beam, p_beam);

    // Small smearing in px and py for emulating "diffractive" pseudorapidity-like proton
    momentum_p1 = CLHEP::HepLorentzVector(px_1, py_1, +p_beam * (1. - xi1), p_beam * (1. - xi1));
    momentum_p2 = CLHEP::HepLorentzVector(px_2, py_2, -p_beam * (1. - xi2), p_beam * (1. - xi2));

    if (verbosity)
    {
      printf("  p_X_z = %.1f\n", momentum_X.z());
      printf("  p_H_z = %.1f\n", momentum_H.z());

      const CLHEP::HepLorentzVector m_tot = momentum_p1 + momentum_p2 + momentum_X + momentum_H;
      printf(" 4-momentum of p + p + X + H: (%.1f, %.1f, %.1f | %.1f)\n", m_tot.x(), m_tot.y(), m_tot.z(), m_tot.t());
    }

    if (momentum_H.perp() > p_T_H_min)
    {
      generationOK = true;
      break;
    }
  }

  if (!generationOK){
    edm::LogWarning("PPSToyMcLHE_XH") << "Failed to generate event.";
    return;
  }

  // Fill in the HepMC record
  unsigned int barcode = 0;

  // Status codes
  const int statusFinal = 1;
  const int statusDecayed = 2;

  int status_X = (decayX) ? statusDecayed : statusFinal;
  int status_H = (decayHToBottoms) ? statusDecayed : statusFinal;

  HepMC::GenParticle* particle_H = new HepMC::GenParticle(momentum_H, particleId_H, status_H);
  particle_H->suggest_barcode(++barcode);
  vtx->add_particle_out(particle_H);

  HepMC::GenParticle* particle_X = new HepMC::GenParticle(momentum_X, particleId_X, status_X);
  particle_X->suggest_barcode(++barcode);
  vtx->add_particle_out(particle_X);

  HepMC::GenParticle* particle_p1 = new HepMC::GenParticle(momentum_p1, particleId_p, statusFinal);
  particle_p1->suggest_barcode(++barcode);
  vtx->add_particle_out(particle_p1);

  HepMC::GenParticle* particle_p2 = new HepMC::GenParticle(momentum_p2, particleId_p, statusFinal);
  particle_p2->suggest_barcode(++barcode);
  vtx->add_particle_out(particle_p2);

  fillEventLHE(hepeup, momentum_p1i, particleId_p, -1, -1);
  fillEventLHE(hepeup, momentum_p2i, particleId_p, -1, -1);
  fillEventLHE(hepeup, momentum_p1, particleId_p, 1, 1);
  fillEventLHE(hepeup, momentum_p2, particleId_p, 1, 1);

  // X and its decays...
  fillEventLHE(hepeup, momentum_S, 51, 2, 0); // DM (S) -> H + X
  fillEventLHE(hepeup, momentum_X, particleId_X, 2, 0);
  if (decayX)
  {
    // generate decay angles in X's rest frame;
    const double cos_theta_d = 2. * CLHEP::RandFlat::shoot(engine) - 1.;
    const double sin_theta_d = sqrt(1. - cos_theta_d * cos_theta_d);
    const double phi_d = CLHEP::RandFlat::shoot(engine) * 2. * M_PI;

    // 4-momentum and energy in X's rest frame
    const double M2 = m_X*m_X - m_f*m_f - m_f*m_f;
    const double p_f = sqrt(M2*M2 - 4. * m_f*m_f * m_f*m_f) / 2. / m_X;
    const double E_f1 = sqrt(p_f*p_f + m_f*m_f);
    const double E_f2 = sqrt(p_f*p_f + m_f*m_f);

    // 4-momenta in X's rest frame
    CLHEP::HepLorentzVector momentum_f1(
	p_f * sin_theta_d * cos(phi_d),
	p_f * sin_theta_d * sin(phi_d),
	p_f * cos_theta_d,
	E_f1
	);

    CLHEP::HepLorentzVector momentum_f2(
	-p_f * sin_theta_d * cos(phi_d),
	-p_f * sin_theta_d * sin(phi_d),
	-p_f * cos_theta_d,
	E_f2
	);

    // apply boost
    double beta = momentum_X.rho() / momentum_X.t();
    CLHEP::Hep3Vector betaVector(momentum_X.x(), momentum_X.y(), momentum_X.z());
    betaVector *= beta / betaVector.mag();
    momentum_f1 = CLHEP::boostOf(momentum_f1, betaVector);
    momentum_f2 = CLHEP::boostOf(momentum_f2, betaVector);

    if (verbosity)
    {
      const CLHEP::HepLorentzVector m_tot = momentum_p1 + momentum_p2 + momentum_H + momentum_f1 + momentum_f2;
      printf("  four-momentum of p + p + H + f1 + f2: (%.5f, %.5f, %.5f | %.5f)\n", m_tot.x(), m_tot.y(), m_tot.z(), m_tot.t());
    }

    // add particles to vertex
    HepMC::GenParticle* particle_f1 = new HepMC::GenParticle(momentum_f1, particleId_f1, statusFinal);
    particle_f1->suggest_barcode(++barcode);
    vtx->add_particle_out(particle_f1);

    HepMC::GenParticle* particle_f2 = new HepMC::GenParticle(momentum_f2, particleId_f2, statusFinal);
    particle_f2->suggest_barcode(++barcode);
    vtx->add_particle_out(particle_f2);

    if(isnan(M2)||isnan(p_f) || isnan(E_f1) ||isnan(E_f2)||isnan(cos_theta_d)||isnan(sin_theta_d)||isnan(phi_d)||isnan(beta)){
      edm::LogWarning("PPSToyMcLHE_XH") << "--NAN--. Skipping the event!";
      return;
    }else{
      fillEventLHE(hepeup, momentum_f1, particleId_f1, 1, 0);
      fillEventLHE(hepeup, momentum_f2, particleId_f2, 1, 0);
    }

  }

  // H and its decays...
  fillEventLHE(hepeup, momentum_H, particleId_H, 2, 0);

  if (decayHToBottoms)
  {
    double m_l = 0.;
    signed int particleId_l_mi=0, particleId_l_pl=0;

    if (decayHToBottoms) m_l = m_b, particleId_l_mi = particleId_b_mi, particleId_l_pl = particleId_b_pl;

    // generate decay angles in H's rest frame;
    const double cos_theta_d = 2. * CLHEP::RandFlat::shoot(engine) - 1.;
    const double sin_theta_d = sqrt(1. - cos_theta_d * cos_theta_d);
    const double phi_d = CLHEP::RandFlat::shoot(engine) * 2. * M_PI;

    // lepton momentum and energy in H's rest frame
    const double E_l = m_H / 2.;
    const double p_l = sqrt(E_l*E_l - m_l*m_l);

    // lepton four-momenta in H's rest frame
    CLHEP::HepLorentzVector momentum_l_mi(
	p_l * sin_theta_d * cos(phi_d),
	p_l * sin_theta_d * sin(phi_d),
	p_l * cos_theta_d,
	E_l
	);

    CLHEP::HepLorentzVector momentum_l_pl(
	-p_l * sin_theta_d * cos(phi_d),
	-p_l * sin_theta_d * sin(phi_d),
	-p_l * cos_theta_d,
	E_l
	);

    // apply boost
    double beta = momentum_H.rho() / momentum_H.t();
    CLHEP::Hep3Vector betaVector(momentum_H.x(), momentum_H.y(), momentum_H.z());
    betaVector *= beta / betaVector.mag();
    momentum_l_mi = CLHEP::boostOf(momentum_l_mi, betaVector);
    momentum_l_pl = CLHEP::boostOf(momentum_l_pl, betaVector);

    if (verbosity)
    {
      const CLHEP::HepLorentzVector m_tot = momentum_p1 + momentum_p2 + momentum_X + momentum_l_mi + momentum_l_pl;
      printf("  four-momentum of p + p + X + q + qbar: (%.1f, %.1f, %.1f | %.1f)\n", m_tot.x(), m_tot.y(), m_tot.z(), m_tot.t());
    }

    // add particles to vertex
    HepMC::GenParticle* particle_l_mi = new HepMC::GenParticle(momentum_l_mi, particleId_l_mi, statusDecayed);
    particle_l_mi->suggest_barcode(++barcode);
    vtx->add_particle_out(particle_l_mi);

    HepMC::GenParticle* particle_l_pl = new HepMC::GenParticle(momentum_l_pl, particleId_l_pl, statusDecayed);
    particle_l_pl->suggest_barcode(++barcode);
    vtx->add_particle_out(particle_l_pl);

    if(isnan(E_l) ||isnan(p_l)||isnan(cos_theta_d)||isnan(sin_theta_d)||isnan(phi_d)||isnan(beta)){
      edm::LogWarning("PPSToyMcLHE_XH") << "--NAN--. Skipping the event!";
      return;
    }else{
      fillEventLHE(hepeup, momentum_l_mi, particleId_l_mi, 1, 0);
      fillEventLHE(hepeup, momentum_l_pl, particleId_l_pl, 1, 0);
    }

  }

  // save output
  std::unique_ptr<HepMCProduct> output(new HepMCProduct()) ;
  output->addHepMCData(fEvt);
  e.put(std::move(output), "unsmeared");

  //double originalXWGTUP_ = 1.;
  double originalXWGTUP_ = 2.8167496e-08;
  std::unique_ptr<LHEEventProduct> LHEevent_(new LHEEventProduct(hepeup, originalXWGTUP_));

  if(writeLHE){
    std::copy(LHEevent_->begin(), LHEevent_->end(), std::ostream_iterator<std::string>(file));
  }

  e.put(std::move(LHEevent_), "source");
}

//----------------------------------------------------------------------------------------------------
LHERunInfoProduct::Header PPSToyMcLHE_XH::SLHA() {

  LHERunInfoProduct::Header exportLHE("slha");
  exportLHE.addLine("\n######################################################################\n");
  exportLHE.addLine("## PARAM_CARD AUTOMATICALY GENERATED BY MG5 FOLLOWING UFO MODEL   ####\n");
  exportLHE.addLine("######################################################################\n");
  exportLHE.addLine("##                                                                  ##\n");
  exportLHE.addLine("##  Width set on Auto will be computed following the information    ##\n");
  exportLHE.addLine("##        present in the decay.py files of the model.               ##\n");
  exportLHE.addLine("##        See  arXiv:1402.1178 for more details.                    ##\n");
  exportLHE.addLine("##                                                                  ##\n");
  exportLHE.addLine("######################################################################\n");
  exportLHE.addLine("\n");
  exportLHE.addLine("###################################\n");
  exportLHE.addLine("## INFORMATION FOR MASS\n");
  exportLHE.addLine("###################################\n");
  exportLHE.addLine("Block mass \n");
  exportLHE.addLine("    6 1.730000e+02 # MT \n");
  exportLHE.addLine("   15 1.777000e+00 # MTA \n");
  exportLHE.addLine("   23 9.118800e+01 # MZ \n");
  exportLHE.addLine("   25 1.250000e+02 # MH \n");
  exportLHE.addLine("## Dependent parameters, given by model restrictions.\n");
  exportLHE.addLine("## Those values should be edited following the \n");
  exportLHE.addLine("## analytical expression. MG5 ignores those values \n");
  exportLHE.addLine("## but they are important for interfacing the output of MG5\n");
  exportLHE.addLine("## to external program such as Pythia.\n");
  exportLHE.addLine("  1 0.000000 # d : 0.0 \n");
  exportLHE.addLine("  2 0.000000 # u : 0.0 \n");
  exportLHE.addLine("  3 0.000000 # s : 0.0 \n");
  exportLHE.addLine("  4 0.000000 # c : 0.0 \n");
  exportLHE.addLine("  5 0.000000 # b : 0.0 \n");
  exportLHE.addLine("  11 0.000000 # e- : 0.0 \n");
  exportLHE.addLine("  12 0.000000 # ve : 0.0 \n");
  exportLHE.addLine("  13 0.000000 # mu- : 0.0 \n");
  exportLHE.addLine("  14 0.000000 # vm : 0.0 \n");
  exportLHE.addLine("  16 0.000000 # vt : 0.0 \n");
  exportLHE.addLine("  21 0.000000 # g : 0.0 \n");
  exportLHE.addLine("  22 0.000000 # a : 0.0 \n");
  exportLHE.addLine(
      "  24 80.419002 # w+ : cmath.sqrt(MZ__exp__2/2. + cmath.sqrt(MZ__exp__4/4. - "
      "(aEW*cmath.pi*MZ__exp__2)/(Gf*sqrt__2))) \n");
  exportLHE.addLine("\n");
  exportLHE.addLine("###################################\n");
  exportLHE.addLine("## INFORMATION FOR SMINPUTS\n");
  exportLHE.addLine("###################################\n");
  exportLHE.addLine("Block sminputs \n");
  exportLHE.addLine("    1 1.325070e+02 # aEWM1 \n");
  exportLHE.addLine("    2 1.166390e-05 # Gf \n");
  exportLHE.addLine("    3 1.180000e-01 # aS \n");
  exportLHE.addLine("\n");
  exportLHE.addLine("###################################\n");
  exportLHE.addLine("## INFORMATION FOR WOLFENSTEIN\n");
  exportLHE.addLine("###################################\n");
  exportLHE.addLine("Block wolfenstein \n");
  exportLHE.addLine("    1 2.253000e-01 # lamWS \n");
  exportLHE.addLine("    2 8.080000e-01 # AWS \n");
  exportLHE.addLine("    3 1.320000e-01 # rhoWS \n");
  exportLHE.addLine("    4 3.410000e-01 # etaWS \n");
  exportLHE.addLine("\n");
  exportLHE.addLine("###################################\n");
  exportLHE.addLine("## INFORMATION FOR YUKAWA\n");
  exportLHE.addLine("###################################\n");
  exportLHE.addLine("Block yukawa \n");
  exportLHE.addLine("    6 1.730000e+02 # ymt \n");
  exportLHE.addLine("   15 1.777000e+00 # ymtau \n");
  exportLHE.addLine("\n");
  exportLHE.addLine("###################################\n");
  exportLHE.addLine("## INFORMATION FOR DECAY\n");
  exportLHE.addLine("###################################\n");
  exportLHE.addLine("DECAY   6 1.491500e+00 # WT \n");
  exportLHE.addLine("DECAY  15 2.270000e-12 # WTau \n");
  exportLHE.addLine("DECAY  23 2.441404e+00 # WZ \n");
  exportLHE.addLine("DECAY  24 2.047600e+00 # WW \n");
  exportLHE.addLine("DECAY  25 6.382339e-03 # WH \n");
  exportLHE.addLine("## Dependent parameters, given by model restrictions.\n");
  exportLHE.addLine("## Those values should be edited following the \n");
  exportLHE.addLine("## analytical expression. MG5 ignores those values \n");
  exportLHE.addLine("## but they are important for interfacing the output of MG5\n");
  exportLHE.addLine("## to external program such as Pythia.\n");
  exportLHE.addLine("DECAY  1 0.000000 # d : 0.0 \n");
  exportLHE.addLine("DECAY  2 0.000000 # u : 0.0 \n");
  exportLHE.addLine("DECAY  3 0.000000 # s : 0.0 \n");
  exportLHE.addLine("DECAY  4 0.000000 # c : 0.0 \n");
  exportLHE.addLine("DECAY  5 0.000000 # b : 0.0 \n");
  exportLHE.addLine("DECAY  11 0.000000 # e- : 0.0 \n");
  exportLHE.addLine("DECAY  12 0.000000 # ve : 0.0 \n");
  exportLHE.addLine("DECAY  13 0.000000 # mu- : 0.0 \n");
  exportLHE.addLine("DECAY  14 0.000000 # vm : 0.0 \n");
  exportLHE.addLine("DECAY  16 0.000000 # vt : 0.0 \n");
  exportLHE.addLine("DECAY  21 0.000000 # g : 0.0 \n");
  exportLHE.addLine("DECAY  22 0.000000 # a : 0.0\n");

  return exportLHE;
}




//----------------------------------------------------------------------------------------------------

PPSToyMcLHE_XH::~PPSToyMcLHE_XH()
{
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(PPSToyMcLHE_XH);
