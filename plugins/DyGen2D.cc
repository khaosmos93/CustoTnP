#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "TH2F.h"
#include "TMath.h"


double Get_Rescale_Weight( double m, std::vector<double> par )
{
  double weight =   par[0]
                  + par[1]*m
                  + par[2]*m*m
                  + par[3]*m*m*m
                  + par[4]*m*m*m*m
                  + par[5]*m*m*m*m*m;

  return weight;
}

class DyGen2D : public edm::EDFilter {
 public:
  explicit DyGen2D(const edm::ParameterSet&);

 private:
 virtual bool filter(edm::Event&, const edm::EventSetup&);

  edm::InputTag src;

  int min_njets;
  int max_njets;
  double min_Y;
  double max_Y;

  double min_pt;
  double max_eta;

  double eventWeight;
  bool useMadgraphWeight;
  double madgraphWeight;
  std::vector<double> rescaleWeights;

  TH2F* Weight_Zmass;
  TH2F* Zpt_Zmass;
  TH2F* Zpz_Zmass;
  TH2F* Zy_Zmass;
  TH2F* Zphi_Zmass;
  TH2F* l_pt_Zmass;
  TH2F* l_eta_Zmass;
  TH2F* l_phi_Zmass;
  TH2F* s_pt_Zmass;
  TH2F* s_eta_Zmass;
  TH2F* s_phi_Zmass;
  TH2F* p_pt_Zmass;
  TH2F* p_eta_Zmass;
  TH2F* p_phi_Zmass;
  TH2F* m_pt_Zmass;
  TH2F* m_eta_Zmass;
  TH2F* m_phi_Zmass;
  TH2F* l_eta_Zy;
  TH2F* s_eta_Zy;
  TH2F* p_eta_Zy;
  TH2F* m_eta_Zy;


  TH2F* Weight_Zmass_;
  TH2F* Zpt_Zmass_;
  TH2F* Zpz_Zmass_;
  TH2F* Zy_Zmass_;
  TH2F* Zphi_Zmass_;
  TH2F* l_pt_Zmass_;
  TH2F* l_eta_Zmass_;
  TH2F* l_phi_Zmass_;
  TH2F* s_pt_Zmass_;
  TH2F* s_eta_Zmass_;
  TH2F* s_phi_Zmass_;
  TH2F* p_pt_Zmass_;
  TH2F* p_eta_Zmass_;
  TH2F* p_phi_Zmass_;
  TH2F* m_pt_Zmass_;
  TH2F* m_eta_Zmass_;
  TH2F* m_phi_Zmass_;
  TH2F* l_eta_Zy_;
  TH2F* s_eta_Zy_;
  TH2F* p_eta_Zy_;
  TH2F* m_eta_Zy_;

};


DyGen2D::DyGen2D(const edm::ParameterSet& cfg)
  : src(cfg.getParameter<edm::InputTag>("src")),
    min_njets(cfg.getParameter<int>("min_njets")),
    max_njets(cfg.getParameter<int>("max_njets")),
    min_Y(cfg.getParameter<double>("min_Y")),
    max_Y(cfg.getParameter<double>("max_Y")),
    min_pt(cfg.getParameter<double>("min_pt")),
    max_eta(cfg.getParameter<double>("max_eta")),
    eventWeight(1.0),
    useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),
    madgraphWeight(1.0),
    rescaleWeights(cfg.getParameter<std::vector<double>>("rescaleWeights"))
{
  consumes<reco::GenParticleCollection>(src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));
  mayConsume<LHEEventProduct>(edm::InputTag("externalLHEProducer"));

  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  Weight_Zmass  = fs->make<TH2F>("Weight_Zmass", "", 100, 0, 10000, 4, -2, 2);
  Zpt_Zmass     = fs->make<TH2F>("Zpt_Zmass", "",    100, 0, 10000, 1000, 0, 10000);
  Zpz_Zmass     = fs->make<TH2F>("Zpz_Zmass", "",    100, 0, 10000, 1000, 0, 10000);
  Zy_Zmass      = fs->make<TH2F>("Zy_Zmass", "",     100, 0, 10000, 96, -4.8, 4.8);
  Zphi_Zmass    = fs->make<TH2F>("Zphi_Zmass", "",   100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  l_pt_Zmass    = fs->make<TH2F>("l_pt_Zmass", "",   100, 0, 10000, 1000, 0, 10000);
  l_eta_Zmass   = fs->make<TH2F>("l_eta_Zmass", "",  100, 0, 10000, 96, -4.8, 4.8);
  l_phi_Zmass   = fs->make<TH2F>("l_phi_Zmass", "",  100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  s_pt_Zmass    = fs->make<TH2F>("s_pt_Zmass", "",   100, 0, 10000, 1000, 0, 10000);
  s_eta_Zmass   = fs->make<TH2F>("s_eta_Zmass", "",  100, 0, 10000, 96, -4.8, 4.8);
  s_phi_Zmass   = fs->make<TH2F>("s_phi_Zmass", "",  100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  p_pt_Zmass    = fs->make<TH2F>("p_pt_Zmass", "",   100, 0, 10000, 1000, 0, 10000);
  p_eta_Zmass   = fs->make<TH2F>("p_eta_Zmass", "",  100, 0, 10000, 96, -4.8, 4.8);
  p_phi_Zmass   = fs->make<TH2F>("p_phi_Zmass", "",  100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  m_pt_Zmass    = fs->make<TH2F>("m_pt_Zmass", "",   100, 0, 10000, 1000, 0, 10000);
  m_eta_Zmass   = fs->make<TH2F>("m_eta_Zmass", "",  100, 0, 10000, 96, -4.8, 4.8);
  m_phi_Zmass   = fs->make<TH2F>("m_phi_Zmass", "",  100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  l_eta_Zy      = fs->make<TH2F>("l_eta_Zy", "",     96, -4.8, 4.8, 96, -4.8, 4.8);
  s_eta_Zy      = fs->make<TH2F>("s_eta_Zy", "",     96, -4.8, 4.8, 96, -4.8, 4.8);
  p_eta_Zy      = fs->make<TH2F>("p_eta_Zy", "",     96, -4.8, 4.8, 96, -4.8, 4.8);
  m_eta_Zy      = fs->make<TH2F>("m_eta_Zy", "",     96, -4.8, 4.8, 96, -4.8, 4.8);

  Weight_Zmass_ = fs->make<TH2F>("Weight_Zmass_", "", 100, 0, 10000, 4, -2, 2);
  Zpt_Zmass_    = fs->make<TH2F>("Zpt_Zmass_", "",    100, 0, 10000, 1000, 0, 10000);
  Zpz_Zmass_    = fs->make<TH2F>("Zpz_Zmass_", "",    100, 0, 10000, 1000, 0, 10000);
  Zy_Zmass_     = fs->make<TH2F>("Zy_Zmass_", "",     100, 0, 10000, 96, -4.8, 4.8);
  Zphi_Zmass_   = fs->make<TH2F>("Zphi_Zmass_", "",   100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  l_pt_Zmass_   = fs->make<TH2F>("l_pt_Zmass_", "",   100, 0, 10000, 1000, 0, 10000);
  l_eta_Zmass_  = fs->make<TH2F>("l_eta_Zmass_", "",  100, 0, 10000, 96, -4.8, 4.8);
  l_phi_Zmass_  = fs->make<TH2F>("l_phi_Zmass_", "",  100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  s_pt_Zmass_   = fs->make<TH2F>("s_pt_Zmass_", "",   100, 0, 10000, 1000, 0, 10000);
  s_eta_Zmass_  = fs->make<TH2F>("s_eta_Zmass_", "",  100, 0, 10000, 96, -4.8, 4.8);
  s_phi_Zmass_  = fs->make<TH2F>("s_phi_Zmass_", "",  100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  p_pt_Zmass_   = fs->make<TH2F>("p_pt_Zmass_", "",   100, 0, 10000, 1000, 0, 10000);
  p_eta_Zmass_  = fs->make<TH2F>("p_eta_Zmass_", "",  100, 0, 10000, 96, -4.8, 4.8);
  p_phi_Zmass_  = fs->make<TH2F>("p_phi_Zmass_", "",  100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  m_pt_Zmass_   = fs->make<TH2F>("m_pt_Zmass_", "",   100, 0, 10000, 1000, 0, 10000);
  m_eta_Zmass_  = fs->make<TH2F>("m_eta_Zmass_", "",  100, 0, 10000, 96, -4.8, 4.8);
  m_phi_Zmass_  = fs->make<TH2F>("m_phi_Zmass_", "",  100, 0, 10000, 50, -TMath::Pi(), TMath::Pi());
  l_eta_Zy_     = fs->make<TH2F>("l_eta_Zy_", "",     96, -4.8, 4.8, 96, -4.8, 4.8);
  s_eta_Zy_     = fs->make<TH2F>("s_eta_Zy_", "",     96, -4.8, 4.8, 96, -4.8, 4.8);
  p_eta_Zy_     = fs->make<TH2F>("p_eta_Zy_", "",     96, -4.8, 4.8, 96, -4.8, 4.8);
  m_eta_Zy_     = fs->make<TH2F>("m_eta_Zy_", "",     96, -4.8, 4.8, 96, -4.8, 4.8);
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


  //-- Get # jets using LHE info
  edm::Handle<LHEEventProduct> LHEInfo;
  event.getByLabel(edm::InputTag("externalLHEProducer"), LHEInfo);

  const lhef::HEPEUP& lheEvent = LHEInfo->hepeup();
  std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;

  int njets    = 0;
  int nleptons = 0;
  for( size_t idxParticle = 0; idxParticle < lheParticles.size(); ++idxParticle ) {
    int id     = lheEvent.IDUP[idxParticle];
    int status = lheEvent.ISTUP[idxParticle];

    if( status != 1)
      continue;

    if( 0 < fabs(id) && fabs(id) < 7 )
      njets++;

    if( 10 < fabs(id) && fabs(id) < 17 )
      nleptons++;
  }

  // std::cout << "njets:    " << njets << std::endl;
  // std::cout << "nleptons: " << nleptons << std::endl;

  if( nleptons != 2 ) {
    edm::LogError("DyGen2D") << "nleptons != 2 in LHE level -> return false";
    return false;
  }

  bool is_njets = njets >= min_njets && njets < max_njets;



  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByLabel(src, genParticles);

  reco::GenParticleCollection::const_iterator genp = genParticles->begin();

  bool isFind1 = false;
  bool isFind2 = false;

  const reco::Candidate *mu1 = 0;
  const reco::Candidate *mu2 = 0;
  reco::Particle::LorentzVector Z;


  bool isFind1_ = false;
  bool isFind2_ = false;

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


  if( isFind1 && isFind2 && isFind1_ && isFind2_ ) {

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

    float p_pt   = -999;
    float p_eta  = -999;
    float p_phi  = -999;
    float m_pt   = -999;
    float m_eta  = -999;
    float m_phi  = -999;

    float p_pt_  = -999;
    float p_eta_ = -999;
    float p_phi_ = -999;
    float m_pt_  = -999;
    float m_eta_ = -999;
    float m_phi_ = -999;

    if( mu1->charge() > mu2->charge() ) {
      p_pt   = mu1->pt();
      p_eta  = mu1->eta();
      p_phi  = mu1->phi();
      m_pt   = mu2->pt();
      m_eta  = mu2->eta();
      m_phi  = mu2->phi();

      p_pt_  = mu1_->pt();
      p_eta_ = mu1_->eta();
      p_phi_ = mu1_->phi();
      m_pt_  = mu2_->pt();
      m_eta_ = mu2_->eta();
      m_phi_ = mu2_->phi();
    }
    else {
      p_pt   = mu2->pt();
      p_eta  = mu2->eta();
      p_phi  = mu2->phi();
      m_pt   = mu1->pt();
      m_eta  = mu1->eta();
      m_phi  = mu1->phi();

      p_pt_  = mu2_->pt();
      p_eta_ = mu2_->eta();
      p_phi_ = mu2_->phi();
      m_pt_  = mu1_->pt();
      m_eta_ = mu1_->eta();
      m_phi_ = mu1_->phi();
    }

    bool is_rapidity = ( fabs(Z.Rapidity()) >= min_Y && fabs(Z.Rapidity()) < max_Y );
    bool is_acc      = ( l_pt_ > min_pt && s_pt_ > min_pt && fabs(l_eta_) < max_eta && fabs(s_eta_) < max_eta );

    bool fill_histo = is_njets && is_rapidity && is_acc;

    //-- Total weight!!!
    double scaleWeight = 1.0;
    if( rescaleWeights.size() == 6 ) {
      scaleWeight = Get_Rescale_Weight( (double)Z_.mass(), rescaleWeights );
    }

    double totalWeight = madgraphWeight * scaleWeight;

    // std::cout << std::endl;
    // std::cout << "scaleWeight: " << scaleWeight << std::endl;
    // std::cout << "totalWeight: " << totalWeight << std::endl;

    if( fill_histo ) {
      Weight_Zmass->Fill( Z.mass(), totalWeight );
      Zpt_Zmass->Fill(    Z.mass(), Z.pt(),  totalWeight );
      Zpz_Zmass->Fill(    Z.mass(), Z.pz(),  totalWeight );
      Zy_Zmass->Fill(     Z.mass(), Z.Rapidity(), totalWeight );
      Zphi_Zmass->Fill(   Z.mass(), Z.phi(), totalWeight );
      l_pt_Zmass->Fill(   Z.mass(), l_pt,    totalWeight );
      l_eta_Zmass->Fill(  Z.mass(), l_eta,   totalWeight );
      l_phi_Zmass->Fill(  Z.mass(), l_phi,   totalWeight );
      s_pt_Zmass->Fill(   Z.mass(), s_pt,    totalWeight );
      s_eta_Zmass->Fill(  Z.mass(), s_eta,   totalWeight );
      s_phi_Zmass->Fill(  Z.mass(), s_phi,   totalWeight );
      p_pt_Zmass->Fill(   Z.mass(), p_pt,    totalWeight );
      p_eta_Zmass->Fill(  Z.mass(), p_eta,   totalWeight );
      p_phi_Zmass->Fill(  Z.mass(), p_phi,   totalWeight );
      m_pt_Zmass->Fill(   Z.mass(), m_pt,    totalWeight );
      m_eta_Zmass->Fill(  Z.mass(), m_eta,   totalWeight );
      m_phi_Zmass->Fill(  Z.mass(), m_phi,   totalWeight );
      l_eta_Zy->Fill(     Z.Rapidity(), l_eta, totalWeight );
      s_eta_Zy->Fill(     Z.Rapidity(), s_eta, totalWeight );
      p_eta_Zy->Fill(     Z.Rapidity(), p_eta, totalWeight );
      m_eta_Zy->Fill(     Z.Rapidity(), m_eta, totalWeight );

      Weight_Zmass_->Fill( Z_.mass(), totalWeight );
      Zpt_Zmass_->Fill(    Z_.mass(), Z_.pt(),  totalWeight );
      Zpz_Zmass_->Fill(    Z_.mass(), Z_.pz(),  totalWeight );
      Zy_Zmass_->Fill(     Z_.mass(), Z_.Rapidity(), totalWeight );
      Zphi_Zmass_->Fill(   Z_.mass(), Z_.phi(), totalWeight );
      l_pt_Zmass_->Fill(   Z_.mass(), l_pt_,    totalWeight );
      l_eta_Zmass_->Fill(  Z_.mass(), l_eta_,   totalWeight );
      l_phi_Zmass_->Fill(  Z_.mass(), l_phi_,   totalWeight );
      s_pt_Zmass_->Fill(   Z_.mass(), s_pt_,    totalWeight );
      s_eta_Zmass_->Fill(  Z_.mass(), s_eta_,   totalWeight );
      s_phi_Zmass_->Fill(  Z_.mass(), s_phi_,   totalWeight );
      p_pt_Zmass_->Fill(   Z_.mass(), p_pt_,    totalWeight );
      p_eta_Zmass_->Fill(  Z_.mass(), p_eta_,   totalWeight );
      p_phi_Zmass_->Fill(  Z_.mass(), p_phi_,   totalWeight );
      m_pt_Zmass_->Fill(   Z_.mass(), m_pt_,    totalWeight );
      m_eta_Zmass_->Fill(  Z_.mass(), m_eta_,   totalWeight );
      m_phi_Zmass_->Fill(  Z_.mass(), m_phi_,   totalWeight );
      l_eta_Zy_->Fill(     Z_.Rapidity(), l_eta_, totalWeight );
      s_eta_Zy_->Fill(     Z_.Rapidity(), s_eta_, totalWeight );
      p_eta_Zy_->Fill(     Z_.Rapidity(), p_eta_, totalWeight );
      m_eta_Zy_->Fill(     Z_.Rapidity(), m_eta_, totalWeight );
    }
  }

  return
    isFind1 && isFind2 && isFind1_ && isFind2_;
}


DEFINE_FWK_MODULE(DyGen2D);
