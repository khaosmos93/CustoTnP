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
  TH2F* Zy_Zmass;
  TH2F* Zphi_Zmass;
  TH2F* l_pt_Zmass;
  TH2F* l_eta_Zmass;
  TH2F* l_phi_Zmass;
  TH2F* s_pt_Zmass;
  TH2F* s_eta_Zmass;
  TH2F* s_phi_Zmass;

  TH2F* Weight_Zmass_;
  TH2F* Zpt_Zmass_;
  TH2F* Zeta_Zmass_;
  TH2F* Zy_Zmass_;
  TH2F* Zphi_Zmass_;
  TH2F* l_pt_Zmass_;
  TH2F* l_eta_Zmass_;
  TH2F* l_phi_Zmass_;
  TH2F* s_pt_Zmass_;
  TH2F* s_eta_Zmass_;
  TH2F* s_phi_Zmass_;

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

  Weight_Zmass  = fs->make<TH2F>("Weight_Zmass", "", 100, 0, 10000, 4, -2, 2);
  Zpt_Zmass     = fs->make<TH2F>("Zpt_Zmass", "",    100, 0, 10000, 10000, 0, 10000);
  Zeta_Zmass    = fs->make<TH2F>("Zeta_Zmass", "",   100, 0, 10000, 200, -10, 10);
  Zy_Zmass      = fs->make<TH2F>("Zy_Zmass", "",     100, 0, 10000, 96, -4.8, 4.8);
  Zphi_Zmass    = fs->make<TH2F>("Zphi_Zmass", "",   100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  l_pt_Zmass    = fs->make<TH2F>("l_pt_Zmass", "",   100, 0, 10000, 10000, 0, 10000);
  l_eta_Zmass   = fs->make<TH2F>("l_eta_Zmass", "",  100, 0, 10000, 96, -4.8, 4.8);
  l_phi_Zmass   = fs->make<TH2F>("l_phi_Zmass", "",  100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  s_pt_Zmass    = fs->make<TH2F>("s_pt_Zmass", "",   100, 0, 10000, 10000, 0, 10000);
  s_eta_Zmass   = fs->make<TH2F>("s_eta_Zmass", "",  100, 0, 10000, 96, -4.8, 4.8);
  s_phi_Zmass   = fs->make<TH2F>("s_phi_Zmass", "",  100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());

  Weight_Zmass_ = fs->make<TH2F>("Weight_Zmass_", "", 100, 0, 10000, 4, -2, 2);
  Zpt_Zmass_    = fs->make<TH2F>("Zpt_Zmass_", "",    100, 0, 10000, 10000, 0, 10000);
  Zeta_Zmass_   = fs->make<TH2F>("Zeta_Zmass_", "",   100, 0, 10000, 200, -10, 10);
  Zy_Zmass_     = fs->make<TH2F>("Zy_Zmass_", "",     100, 0, 10000, 96, -4.8, 4.8);
  Zphi_Zmass_   = fs->make<TH2F>("Zphi_Zmass_", "",   100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  l_pt_Zmass_   = fs->make<TH2F>("l_pt_Zmass_", "",   100, 0, 10000, 10000, 0, 10000);
  l_eta_Zmass_  = fs->make<TH2F>("l_eta_Zmass_", "",  100, 0, 10000, 96, -4.8, 4.8);
  l_phi_Zmass_  = fs->make<TH2F>("l_phi_Zmass_", "",  100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  s_pt_Zmass_   = fs->make<TH2F>("s_pt_Zmass_", "",   100, 0, 10000, 10000, 0, 10000);
  s_eta_Zmass_  = fs->make<TH2F>("s_eta_Zmass_", "",  100, 0, 10000, 96, -4.8, 4.8);
  s_phi_Zmass_  = fs->make<TH2F>("s_phi_Zmass_", "",  100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
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
  const reco::Candidate *mu1 = 0;
  const reco::Candidate *mu2 = 0;
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
  const reco::Candidate *mu1_ = 0;
  const reco::Candidate *mu2_ = 0;
  reco::Particle::LorentzVector Z_;


  for (; genp != genParticles->end(); genp++) {

    if( (genp->pdgId() == 13) && (genp->isHardProcess()) ){
      isFind1 = true;
      mu1 = &*genp;
    }

    if( (genp->pdgId() == -13) && (genp->isHardProcess()) ){
      isFind2 = true;
      mu2 = &*genp;
    }

    if( (genp->pdgId() == 13) && (genp->status() == 1) ){

      const reco::Candidate* m = genp->mother();
      bool ok = false;
      while (m) {
        if( m->pdgId() == 32 || m->pdgId() == 23 || m->pdgId() == 39 || m->pdgId() == 5000039 ) {
          ok = true;
          break;
        }
        m = m->mother();
      }

      if(ok) {
        isFind1_ = true;
        mu1_ = &*genp;
      }
    }

    if( (genp->pdgId() == -13) && (genp->status() == 1) ){

      const reco::Candidate* m = genp->mother();
      bool ok = false;
      while (m) {
        if( m->pdgId() == 32 || m->pdgId() == 23 || m->pdgId() == 39 || m->pdgId() == 5000039 ) {
          ok = true;
          break;
        }
        m = m->mother();
      }

      if(ok) {
        isFind2_ = true;
        mu2_ = &*genp;
      }
    }

  }


  if( isFind1 && isFind2 && isFind1_ && isFind2_ ){

    Z  = mu1->p4() + mu2->p4();
    Z_ = mu1_->p4() + mu2_->p4();

    float l_pt   = -999;
    float l_eta  = -999;
    float l_phi  = -999;
    float s_pt   = -999;
    float s_eta  = -999;
    float s_phi  = -999;

    float l_pt_  = -999;
    float l_eta_ = -999;
    float l_phi_ = -999;
    float s_pt_  = -999;
    float s_eta_ = -999;
    float s_phi_ = -999;

    if( mu1_->pt() > mu2_->pt() ) {
      l_pt   = mu1->pt();
      l_eta  = mu1->eta();
      l_phi  = mu1->phi();
      s_pt   = mu2->pt();
      s_eta  = mu2->eta();
      s_phi  = mu2->phi();

      l_pt_  = mu1_->pt();
      l_eta_ = mu1_->eta();
      l_phi_ = mu1_->phi();
      s_pt_  = mu2_->pt();
      s_eta_ = mu2_->eta();
      s_phi_ = mu2_->phi();
    }
    else {
      l_pt   = mu2->pt();
      l_eta  = mu2->eta();
      l_phi  = mu2->phi();
      s_pt   = mu1->pt();
      s_eta  = mu1->eta();
      s_phi  = mu1->phi();

      l_pt_  = mu2_->pt();
      l_eta_ = mu2_->eta();
      l_phi_ = mu2_->phi();
      s_pt_  = mu1_->pt();
      s_eta_ = mu1_->eta();
      s_phi_ = mu1_->phi();
    }

    Weight_Zmass->Fill( Z.mass(), madgraphWeight );
    Zpt_Zmass->Fill(    Z.mass(), Z.pt(),  madgraphWeight );
    Zeta_Zmass->Fill(   Z.mass(), Z.eta(), madgraphWeight );
    Zy_Zmass->Fill(     Z.mass(), Z.Rapidity(), madgraphWeight );
    Zphi_Zmass->Fill(   Z.mass(), Z.phi(), madgraphWeight );
    l_pt_Zmass->Fill(   Z.mass(), l_pt,    madgraphWeight );
    l_eta_Zmass->Fill(  Z.mass(), l_eta,   madgraphWeight );
    l_phi_Zmass->Fill(  Z.mass(), l_phi,   madgraphWeight );
    s_pt_Zmass->Fill(   Z.mass(), s_pt,    madgraphWeight );
    s_eta_Zmass->Fill(  Z.mass(), s_eta,   madgraphWeight );
    s_phi_Zmass->Fill(  Z.mass(), s_phi,   madgraphWeight );

    Weight_Zmass_->Fill( Z_.mass(), madgraphWeight );
    Zpt_Zmass_->Fill(    Z_.mass(), Z_.pt(),  madgraphWeight );
    Zeta_Zmass_->Fill(   Z_.mass(), Z_.eta(), madgraphWeight );
    Zy_Zmass_->Fill(     Z_.mass(), Z_.Rapidity(), madgraphWeight );
    Zphi_Zmass_->Fill(   Z_.mass(), Z_.phi(), madgraphWeight );
    l_pt_Zmass_->Fill(   Z_.mass(), l_pt_,    madgraphWeight );
    l_eta_Zmass_->Fill(  Z_.mass(), l_eta_,   madgraphWeight );
    l_phi_Zmass_->Fill(  Z_.mass(), l_phi_,   madgraphWeight );
    s_pt_Zmass_->Fill(   Z_.mass(), s_pt_,    madgraphWeight );
    s_eta_Zmass_->Fill(  Z_.mass(), s_eta_,   madgraphWeight );
    s_phi_Zmass_->Fill(  Z_.mass(), s_phi_,   madgraphWeight );

  }

  return
    isFind1 && isFind2 && isFind1_ && isFind2_;
}


DEFINE_FWK_MODULE(DyGen2D);
