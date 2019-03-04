#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TH2F.h"
#include "TMath.h"


class DyGen2D : public edm::EDFilter {
 public:
  explicit DyGen2D(const edm::ParameterSet&);

 private:
 virtual bool filter(edm::Event&, const edm::EventSetup&);

  edm::InputTag src;

  double eventWeight;
  bool useMadgraphWeight;
  double madgraphWeight;

  TH2F* Weight_Zmass;
  TH2F* Zpt_Zmass;
  TH2F* Zeta_Zmass;
  TH2F* Zphi_Zmass;
  TH2F* pt_Zmass;
  TH2F* eta_Zmass;
  TH2F* phi_Zmass;

  TH2F* Weight_Zmass_;
  TH2F* Zpt_Zmass_;
  TH2F* Zeta_Zmass_;
  TH2F* Zphi_Zmass_;
  TH2F* pt_Zmass_;
  TH2F* eta_Zmass_;
  TH2F* phi_Zmass_;

};


DyGen2D::DyGen2D(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),
    eventWeight(1.0),
    useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),
    madgraphWeight(1.0)
{
  consumes<reco::GenParticleCollection>(src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));

  edm::Service<TFileService> fs;

  Weight_Zmass = fs->make<TH2F>("Weight_Zmass", "", 100, 0, 10000, 4, -2, 2);
  Zpt_Zmass    = fs->make<TH2F>("Zpt_Zmass", "",    100, 0, 10000, 10000, 0, 10000);
  Zeta_Zmass   = fs->make<TH2F>("Zeta_Zmass", "",   100, 0, 10000, 96, -4.8, 4.8);
  Zphi_Zmass   = fs->make<TH2F>("Zphi_Zmass", "",   100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  pt_Zmass     = fs->make<TH2F>("pt_Zmass", "",     100, 0, 10000, 10000, 0, 10000);
  eta_Zmass    = fs->make<TH2F>("eta_Zmass", "",    100, 0, 10000, 96, -4.8, 4.8);
  phi_Zmass    = fs->make<TH2F>("phi_Zmass", "",    100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());

  Weight_Zmass_ = fs->make<TH2F>("Weight_Zmass_", "", 100, 0, 10000, 4, -2, 2);
  Zpt_Zmass_    = fs->make<TH2F>("Zpt_Zmass_", "",    100, 0, 10000, 10000, 0, 10000);
  Zeta_Zmass_   = fs->make<TH2F>("Zeta_Zmass_", "",   100, 0, 10000, 96, -4.8, 4.8);
  Zphi_Zmass_   = fs->make<TH2F>("Zphi_Zmass_", "",   100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  pt_Zmass_     = fs->make<TH2F>("pt_Zmass_", "",     100, 0, 10000, 10000, 0, 10000);
  eta_Zmass_    = fs->make<TH2F>("eta_Zmass_", "",    100, 0, 10000, 96, -4.8, 4.8);
  phi_Zmass_    = fs->make<TH2F>("phi_Zmass_", "",    100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
}

bool DyGen2D::filter(edm::Event& event, const edm::EventSetup&) {

  eventWeight = 1.0;
  madgraphWeight = 1.0;

  if (useMadgraphWeight) {
    edm::Handle<GenEventInfoProduct> gen_ev_info;
    event.getByLabel(edm::InputTag("generator"), gen_ev_info);
    if (gen_ev_info.isValid() ){
      eventWeight = gen_ev_info->weight();
      madgraphWeight = ( eventWeight > 0 ) ? 1.0 : -1.0;
    }
  }
  else {
    eventWeight = 1.0;
    madgraphWeight = 1.0;
  }

  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByLabel(src, genParticles);

  reco::GenParticleCollection::const_iterator genp = genParticles->begin();

  bool isFind1 = false;
  bool isFind2 = false;

  // float pt1 = 0;
  // float pt2 = 0;
  // float eta1 =0;
  // float eta2 =0;
  // float phi1 =0;
  // float phi2 =0;
  // float mass1 = 0;
  // float mass2 = 0;
  reco::Particle::LorentzVector mu1, mu2;
  reco::Particle::LorentzVector Z;


  bool isFind1_ = false;
  bool isFind2_ = false;

  // float pt1_ = 0;
  // float pt2_ = 0;
  // float eta1_ =0;
  // float eta2_ =0;
  // float phi1_ =0;
  // float phi2_ =0;
  // float mass1_ = 0;
  // float mass2_ = 0;
  reco::Particle::LorentzVector mu1_, mu2_;
  reco::Particle::LorentzVector Z_;


  for (; genp != genParticles->end(); genp++) {

    if( (genp->pdgId() == 13) && (genp->isHardProcess()) ){
      isFind1 = true;
      mu1 = genp->p4();
      // eta1    = genp->eta();
      // phi1    = genp->phi();
      // pt1     = genp->pt();
      // mass1   = genp->mass();
      // mu1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);  //0.105
    }

    if( (genp->pdgId() == -13) && (genp->isHardProcess()) ){
      isFind2 = true;
      mu2 = genp->p4();
      // eta2    = genp->eta();
      // phi2    = genp->phi();
      // pt2     = genp->pt();
      // mass2   = genp->mass();
      // mu2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);  //0.105
    }

    if( (genp->pdgId() == 13) && (genp->status() == 1) ){
      isFind1_ = true;
      mu1_ = genp->p4();
      // eta1_    = genp->eta();
      // phi1_    = genp->phi();
      // pt1_     = genp->pt();
      // mass1_   = genp->mass();
      // mu1_.SetPtEtaPhiM(pt1_, eta1_, phi1_, mass1_);  //0.105
    }

    if( (genp->pdgId() == -13) && (genp->status() == 1) ){
      isFind2_ = true;
      mu2_ = genp->p4();
      // eta2_    = genp->eta();
      // phi2_    = genp->phi();
      // pt2_     = genp->pt();
      // mass2_   = genp->mass();
      // mu2_.SetPtEtaPhiM(pt2_, eta2_, phi2_, mass2_);  //0.105
    }

  }

  Z = mu1 + mu2;

  Z_ = mu1_ + mu2_;


  if( isFind1 && isFind2 && isFind1_ && isFind2_ ){

    if(Z_.mass() < 100) {
      std::cout << "M: " <<  << std::endl;
      std::cout << "\t p+: " << mu1_ << std::endl;
      std::cout << "\t p-: " << mu2_ << std::endl;
      std::cout << "\t\t p4: " << (mu1_+mu2_) << std::endl;
    }

    float l_pt   = mu1_.pt() > mu2_.pt() ? mu1.pt()  : mu2.pt();
    float l_eta  = mu1_.pt() > mu2_.pt() ? mu1.eta() : mu2.eta();
    float l_phi  = mu1_.pt() > mu2_.pt() ? mu1.phi() : mu2.phi();

    float l_pt_   = mu1_.pt() > mu2_.pt() ? mu1_.pt()  : mu2_.pt();
    float l_eta_  = mu1_.pt() > mu2_.pt() ? mu1_.eta() : mu2_.eta();
    float l_phi_  = mu1_.pt() > mu2_.pt() ? mu1_.phi() : mu2_.phi();

    Weight_Zmass->Fill( Z.mass(), madgraphWeight );
    Zpt_Zmass->Fill(    Z.mass(), Z.pt(),  madgraphWeight );
    Zeta_Zmass->Fill(   Z.mass(), Z.eta(), madgraphWeight );
    Zphi_Zmass->Fill(   Z.mass(), Z.phi(), madgraphWeight );
    pt_Zmass->Fill(     Z.mass(), l_pt,    madgraphWeight );
    eta_Zmass->Fill(    Z.mass(), l_eta,   madgraphWeight );
    phi_Zmass->Fill(    Z.mass(), l_phi,   madgraphWeight );

    Weight_Zmass_->Fill( Z_.mass(), madgraphWeight );
    Zpt_Zmass_->Fill(    Z_.mass(), Z_.pt(),  madgraphWeight );
    Zeta_Zmass_->Fill(   Z_.mass(), Z_.eta(), madgraphWeight );
    Zphi_Zmass_->Fill(   Z_.mass(), Z_.phi(), madgraphWeight );
    pt_Zmass_->Fill(     Z_.mass(), l_pt_,    madgraphWeight );
    eta_Zmass_->Fill(    Z_.mass(), l_eta_,   madgraphWeight );
    phi_Zmass_->Fill(    Z_.mass(), l_phi_,   madgraphWeight );

  }

  return
    isFind1 && isFind2 && isFind1_ && isFind2_;
}


DEFINE_FWK_MODULE(DyGen2D);
