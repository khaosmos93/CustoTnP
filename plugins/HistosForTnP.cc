#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CustoTnP/Analyzer/src/DileptonUtilities.h"
#include "CustoTnP/Analyzer/src/GeneralUtilities.h"
#include "CustoTnP/Analyzer/src/PATUtilities.h"
#include "CustoTnP/Analyzer/src/ToConcrete.h"
#include "CustoTnP/Analyzer/src/TrackUtilities.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//using namespace std;
void CopyVectorToArray( std::vector<double>& vec, double*& arr ) {
  int nBin = (int)vec.size()-1; // -- # bins = # bin edges-1 -- //
  arr = new Double_t[nBin+1]; // -- dynamic allocation -- //
  for(int i=0; i<nBin+1; i++)
    arr[i] = vec[i];
}

class CustoTnPHistosForTnP : public edm::EDAnalyzer {
 public:
  explicit CustoTnPHistosForTnP(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

 private:
  void getBSandPV(const edm::Event&);
  void fillTnPControlHistos(const pat::CompositeCandidate&, const reco::CandidateBaseRef&, const reco::CandidateBaseRef&, float, int, bool );

  typedef std::pair< std::vector<TH1F*>, std::vector<TH1F*> > BinHistos;
  void fillTnPBinHistos( double, double, bool, double, std::vector<double>& , BinHistos&, bool );

  typedef std::vector< TH2F* > BinHistos2D;
  void fillTnPBinHistos2D( double, double, bool, double, BinHistos2D&, bool );

  int calcNShowers(const reco::CandidateBaseRef&);

  edm::InputTag dilepton_src;
  edm::InputTag beamspot_src;
  edm::InputTag vertex_src;
  const bool use_bs_and_pv;
  const reco::BeamSpot* beamspot;
  const reco::Vertex*   vertex;
  int                   nVtx;

  StringCutObjectSelector<pat::Muon> passing_probe_selector;
  StringCutObjectSelector<pat::Muon> comparison_probe_selector;

  double minMass;
  double maxMass;
  double probe_pt_min;
  double passing_probe_dpt_over_pt_max;
  double passing_probe_dz_max;
  double comparison_probe_dpt_over_pt_max;
  double comparison_probe_dz_max;

  bool   _usePrescaleWeight;
  int    _prescaleWeight;
  double _eventWeight;
  bool   _useMadgraphWeight;
  double _madgraphWeight;

  edm::InputTag pileup_src;
  double _pileupWeight;
  std::vector<double> vec_PileUpWeight;

  double _totalWeight;

  bool isAOD;
  bool useBinHistos2D;
  const bool ShutUp;


  // -- Histograms -- //
  TH1F* NBeamSpot;
  TH1F* NVertices;
  TH1F* NDileptons;
  TH1F* WeightMadGraph;

  TH1F *TagPt;
  TH1F *TagEta;
  TH1F *TagPhi;

  /*
    TH1F* TagAbsTkIso;
    TH1F* TagRelTkIso;
    TH1F* TagChi2dof;
    TH1F* TagTrackDXYBS;
    TH1F* TagTrackDZBS;
    TH1F* TagTrackDXYPV;
    TH1F* TagTrackDZPV;
    TH1F* TagNPxHits;
    TH1F* TagNStHits;
    TH1F* TagNTkHits;
    TH1F* TagNMuHits;
    TH1F* TagNHits;
    TH1F* TagNPxLayers;
    TH1F* TagNStLayers;
    TH1F* TagNTkLayers;
  */

  TH1F *ProbePt;
  TH1F *ProbeEta;
  TH1F *ProbePhi;
  TH1F *ProbeNVertices;
  TH2F *ProbeEtaPhi;
  TH2F *ProbeEtaPt;
  TH2F *ProbeEtaDptPt;
  TH2F *ProbeEtaShower;

  TH1F *PassingProbePt;
  TH1F *PassingProbeEta;
  TH1F *PassingProbePhi;
  TH1F *PassingProbeNVertices;
  TH2F *PassingProbeEtaPhi;
  TH2F *PassingProbeEtaPt;
  TH2F *PassingProbeEtaDptPt;
  TH2F *PassingProbeEtaShower;

  TH1F *FailingProbePt;
  TH1F *FailingProbeEta;
  TH1F *FailingProbePhi;
  TH1F *FailingProbeNVertices;
  TH2F *FailingProbeEtaPhi;
  TH2F *FailingProbeEtaPt;
  TH2F *FailingProbeEtaDptPt;
  TH2F *FailingProbeEtaShower;

  TH1F *PairNoPtMass;
  TH1F *PairNoPtPt;
  TH1F *PairNoPtEta;
  TH1F *PairNoPtRap;

  TH1F *PairMass;
  TH1F *PairPt;
  TH1F *PairEta;
  TH1F *PairRap;

  TH1F *PassingPairMass;
  TH1F *PassingPairPt;
  TH1F *PassingPairEta;
  TH1F *PassingPairRap;

  TH1F *FailingPairMass;
  TH1F *FailingPairPt;
  TH1F *FailingPairEta;
  TH1F *FailingPairRap;

  /*
    TH1F* ProbeAbsTkIso;
    TH1F* ProbeRelTkIso;
    TH1F* ProbeChi2dof;
    TH1F* ProbeTrackDXYBS;
    TH1F* ProbeTrackDZBS;
    TH1F* ProbeTrackDXYPV;
    TH1F* ProbeTrackDZPV;
    TH1F* ProbeNPxHits;
    TH1F* ProbeNStHits;
    TH1F* ProbeNTkHits;
    TH1F* ProbeNMuHits;
    TH1F* ProbeNHits;
    TH1F* ProbeNPxLayers;
    TH1F* ProbeNStLayers;
    TH1F* ProbeNTkLayers;
  */

  // -- bin distributions for passing and failing probes -- //
  std::vector<double> vec_PtBins;
  std::vector<double> vec_AbsPBins;
  std::vector<double> vec_EtaBins;
  std::vector<double> vec_PhiBins;
  std::vector<double> vec_VtxBins;
  std::vector<double> vec_ShowerBins;

  TH1F *hEffTemplatePt;
  TH1F *hEffTemplateAbsP;
  TH1F *hEffTemplateEta;
  TH1F *hEffTemplatePhi;
  TH1F *hEffTemplateVtx;
  TH1F *hEffTemplateShower;

  std::pair< std::vector<TH1F*>, std::vector<TH1F*> > make_bin_histos(TString name, std::vector<double>& vec_bins) {
    edm::Service<TFileService> fs;

    int nbins = ( (int)vec_bins.size() ) - 1;

    unsigned nMassBin = (unsigned)(maxMass - minMass);

    std::vector<TH1F*> vec_pass;
    std::vector<TH1F*> vec_fail;

    TString histNameBase = TString::Format("h%s", name.Data());

    for(int i=0; i<nbins; ++i) {
      TString binInfo = TString::Format("%02dbin", i);

      TString histNamePass = histNameBase + "Pass_" + binInfo;
      TH1F* hTempPass = fs->make<TH1F>(histNamePass, "", nMassBin, minMass, maxMass);  hTempPass->Sumw2();
      vec_pass.push_back( hTempPass );

      TString histNameFail = histNameBase + "Fail_" + binInfo;
      TH1F* hTempFail = fs->make<TH1F>(histNameFail, "", nMassBin, minMass, maxMass);  hTempFail->Sumw2();
      vec_fail.push_back( hTempFail );
    }

    return std::make_pair(vec_pass, vec_fail);
  }

  std::vector< TH2F* > make_bin_histos_2D(TString name, std::vector<double>& vec_bins) {
    edm::Service<TFileService> fs;

    std::vector< TH2F* > vec_h2d = {};

    Double_t mass_bins_for_2D[13] = {0, 60, 120, 200, 300, 400, 800, 1400, 2300, 3500, 4500, 6000, 1e6};
    unsigned nMassBin = (unsigned)(12);

    int nbins = ( (int)vec_bins.size() ) - 1;
    double *arr_bins;
    CopyVectorToArray(vec_bins, arr_bins);

    TString histNameBase = TString::Format("h2d%s", name.Data());

    TH2F* h2d_pass = fs->make<TH2F>(histNameBase + "Pass", "", nMassBin, mass_bins_for_2D, nbins, arr_bins );  h2d_pass->Sumw2();
    TH2F* h2d_fail = fs->make<TH2F>(histNameBase + "Fail", "", nMassBin, mass_bins_for_2D, nbins, arr_bins );  h2d_fail->Sumw2();
    TH2F* h2d_sum  = fs->make<TH2F>(histNameBase + "Sum",  "", nMassBin, mass_bins_for_2D, nbins, arr_bins );  h2d_sum->Sumw2();
    TH2F* h2d_sqs  = fs->make<TH2F>(histNameBase + "Sqs",  "", nMassBin, mass_bins_for_2D, nbins, arr_bins );  h2d_sqs->Sumw2();

    vec_h2d.push_back( h2d_pass );
    vec_h2d.push_back( h2d_fail );
    vec_h2d.push_back( h2d_sum );
    vec_h2d.push_back( h2d_sqs );

    return vec_h2d;
  }

  // -- Bin histograms -- //
  BinHistos Pt;
  BinHistos PtB;  // 0.0 < |eta| < 0.9
  BinHistos PtO;  // 0.9 < |eta| < 1.2
  BinHistos PtE;  // 1.2 < |eta| < 2.1
  BinHistos PtF;  // 2.1 < |eta| < 2.4

  BinHistos AbsP;
  BinHistos AbsPB;  // 0.0 < |eta| < 0.9
  BinHistos AbsPO;  // 0.9 < |eta| < 1.2
  BinHistos AbsPE;  // 1.2 < |eta| < 2.1
  BinHistos AbsPF;  // 2.1 < |eta| < 2.4

  BinHistos Eta;
  BinHistos Phi;
  BinHistos Vtx;

  BinHistos Shower;
  BinHistos ShowerE;
  BinHistos ShowerB;

  BinHistos2D Pt2D;
  BinHistos2D PtB2D;  // 0.0 < |eta| < 0.9
  BinHistos2D PtO2D;  // 0.9 < |eta| < 1.2
  BinHistos2D PtE2D;  // 1.2 < |eta| < 2.1
  BinHistos2D PtF2D;  // 2.1 < |eta| < 2.4

  BinHistos2D AbsP2D;
  BinHistos2D AbsPB2D;  // 0.0 < |eta| < 0.9
  BinHistos2D AbsPO2D;  // 0.9 < |eta| < 1.2
  BinHistos2D AbsPE2D;  // 1.2 < |eta| < 2.1
  BinHistos2D AbsPF2D;  // 2.1 < |eta| < 2.4

  BinHistos2D Eta2D;
  BinHistos2D Phi2D;
  BinHistos2D Vtx2D;

  BinHistos2D Shower2D;
  BinHistos2D ShowerE2D;
  BinHistos2D ShowerB2D;

  bool IsRealData;
  int RunNum;
  int LumiBlockNum;
  unsigned long long EventNum;
  double Mass;
  double VertexMass;
  double Probe_Pt;
  double Probe_Eta;
  TTree* comparison_tree;
};

CustoTnPHistosForTnP::CustoTnPHistosForTnP(const edm::ParameterSet& cfg)
  : dilepton_src(cfg.getParameter<edm::InputTag>("dilepton_src")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    use_bs_and_pv(cfg.getParameter<bool>("use_bs_and_pv")),
    beamspot(0),
    vertex(0),
    nVtx(0),

    passing_probe_selector(cfg.getParameter<std::string>("passing_probe_cut")),
    comparison_probe_selector(cfg.getParameter<std::string>("comparison_probe_cut")),

    minMass(cfg.getParameter<double>("minMass")),
    maxMass(cfg.getParameter<double>("maxMass")),

    probe_pt_min(cfg.getParameter<double>("probe_pt_min")),
    passing_probe_dpt_over_pt_max(cfg.getParameter<double>("passing_probe_dpt_over_pt_max")),
    passing_probe_dz_max(cfg.getParameter<double>("passing_probe_dz_max")),
    comparison_probe_dpt_over_pt_max(cfg.getParameter<double>("comparison_probe_dpt_over_pt_max")),
    comparison_probe_dz_max(cfg.getParameter<double>("comparison_probe_dz_max")),

    _usePrescaleWeight(cfg.getUntrackedParameter<bool>("usePrescaleWeight",false)),
    _prescaleWeight(1),
    _useMadgraphWeight(cfg.getParameter<bool>("useMadgraphWeight")),
    _madgraphWeight(1.),

    pileup_src(cfg.getParameter<edm::InputTag>("pileup_src")),
    _pileupWeight(1.),
    vec_PileUpWeight(cfg.getParameter<std::vector<double>>("vec_PileUpWeight")),

    _totalWeight(1.),

    isAOD(cfg.getParameter<bool>("isAOD")),
    useBinHistos2D(cfg.getParameter<bool>("useBinHistos2D")),

    ShutUp(cfg.getParameter<bool>("ShutUp")),

    vec_PtBins(cfg.getParameter<std::vector<double>>("vec_PtBins")),
    vec_AbsPBins(cfg.getParameter<std::vector<double>>("vec_AbsPBins")),
    vec_EtaBins(cfg.getParameter<std::vector<double>>("vec_EtaBins")),
    vec_PhiBins(cfg.getParameter<std::vector<double>>("vec_PhiBins")),
    vec_VtxBins(cfg.getParameter<std::vector<double>>("vec_VtxBins")),
    vec_ShowerBins(cfg.getParameter<std::vector<double>>("vec_ShowerBins")),
{
  consumes<pat::CompositeCandidateCollection>(dilepton_src);
  consumes<reco::BeamSpot>(beamspot_src);
  consumes<reco::VertexCollection>(vertex_src);
  consumes<std::vector<PileupSummaryInfo> >(pileup_src);
  mayConsume<GenEventInfoProduct>(edm::InputTag("generator"));

  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  NBeamSpot = fs->make<TH1F>("NBeamSpot", "# beamspots/event",  2, 0,  2);
  NVertices = fs->make<TH1F>("NVertices", "# vertices/event",  200, 0, 200);

  // Dilepton multiplicity.
  NDileptons = fs->make<TH1F>("NDileptons", "# dileptons/event", 10, 0, 10);

  // Generator weights
  WeightMadGraph = fs->make<TH1F>("WeightMadGraph", "weight per event", 4, -2,2);

  // Tag
  TagPt = fs->make<TH1F>("TagPt", "Tag pT", 10000, 0, 10000);
  TagEta = fs->make<TH1F>("TagEta", "Tag #eta",    96, -4.8, 4.8);
  TagPhi = fs->make<TH1F>("TagPhi", "Tag #phi", 41, -TMath::Pi(), TMath::Pi());

  /*
    TagAbsTkIso = fs->make<TH1F>("TagAbsTkIso", "Tag Iso. (#Delta R < 0.3) #Sigma pT", 1000, 0, 1000);
    TagRelTkIso = fs->make<TH1F>("TagRelTkIso", "Tag Iso. (#Delta R < 0.3) #Sigma pT / tk. pT", 500, 0, 5);
    TagChi2dof = fs->make<TH1F>("TagChi2dof", "Tag #chi^{2}/dof", 500, 0, 50);
    TagTrackDXYBS = fs->make<TH1F>("TagTrackDXYBS", "Tag |dxy wrt BS|", 10000, 0, 2);
    TagTrackDZBS = fs->make<TH1F>("TagTrackDZBS", "Tag |dz wrt BS|", 10000, 0, 20);
    TagTrackDXYPV = fs->make<TH1F>("TagTrackDXYPV", "Tag |dxy wrt PV|", 10000, 0, 2);
    TagTrackDZPV = fs->make<TH1F>("TagTrackDZPV", "Tag |dz wrt PV|", 10000, 0, 20);
    TagNPxHits = fs->make<TH1F>("TagNPxHits", "Tag # pixel hits", 10, 0,  10);
    TagNStHits = fs->make<TH1F>("TagNStHits", "Tag # strip hits", 40, 0, 40);
    TagNTkHits = fs->make<TH1F>("TagNTkHits", "Tag # tracker hits", 50, 0, 50);
    TagNMuHits = fs->make<TH1F>("TagNMuHits", "Tag # muon hits", 60, 0, 60);
    TagNHits = fs->make<TH1F>("TagNHits", "Tag # hits", 80, 0, 80);
    TagNPxLayers = fs->make<TH1F>("TagNPxLayers", "Tag # pixel layers", 10, 0, 10);
    TagNStLayers = fs->make<TH1F>("TagNStLayers", "Tag # strip layers", 20, 0, 20);
    TagNTkLayers = fs->make<TH1F>("TagNTkLayers", "Tag # tracker layers", 30, 0, 30);
  */

  Double_t eta_bins_for_2D[15] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.0, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4};

  // Probe
  ProbePt               = fs->make<TH1F>("ProbePt", "Probe pT", 10000, 0, 10000);
  ProbeEta              = fs->make<TH1F>("ProbeEta", "Probe #eta",    96, -4.8, 4.8);
  ProbePhi              = fs->make<TH1F>("ProbePhi", "Probe #phi", 41, -TMath::Pi(), TMath::Pi());
  ProbeNVertices        = fs->make<TH1F>("ProbeNVertices", "Probe # vertices/event",  200, 0, 200);
  ProbeEtaPhi           = fs->make<TH2F>("ProbeEtaPhi", "Probe #eta #phi",    14, eta_bins_for_2D, 41, -TMath::Pi(), TMath::Pi());
  ProbeEtaPt            = fs->make<TH2F>("ProbeEtaPt", "Probe #eta pT",    14, eta_bins_for_2D, 10000, 0, 10000);
  ProbeEtaDptPt         = fs->make<TH2F>("ProbeEtaDptPt", "Probe #eta dpt/pT",    14, eta_bins_for_2D, 1000, 0, 2);
  if(isAOD)
    ProbeEtaShower        = fs->make<TH2F>("ProbeEtaShower", "Probe #eta shower",    14, eta_bins_for_2D, 7, -1.5, 5.5);

  PassingProbePt        = fs->make<TH1F>("PassingProbePt", "PassingProbe pT", 10000, 0, 10000);
  PassingProbeEta       = fs->make<TH1F>("PassingProbeEta", "PassingProbe #eta",    96, -4.8, 4.8);
  PassingProbePhi       = fs->make<TH1F>("PassingProbePhi", "PassingProbe #phi", 41, -TMath::Pi(), TMath::Pi());
  PassingProbeNVertices = fs->make<TH1F>("PassingProbeNVertices", "PassingProbe # vertices/event",  200, 0, 200);
  PassingProbeEtaPhi    = fs->make<TH2F>("PassingProbeEtaPhi", "PassingProbe #eta #phi",    14, eta_bins_for_2D, 41, -TMath::Pi(), TMath::Pi());
  PassingProbeEtaPt     = fs->make<TH2F>("PassingProbeEtaPt", "PassingProbe #eta pT",    14, eta_bins_for_2D, 10000, 0, 10000);
  PassingProbeEtaDptPt  = fs->make<TH2F>("PassingProbeEtaDptPt", "PassingProbe #eta dpt/pT",    14, eta_bins_for_2D, 1000, 0, 2);
  if(isAOD)
    PassingProbeEtaShower = fs->make<TH2F>("PassingProbeEtaShower", "PassingProbe #eta shower",    14, eta_bins_for_2D, 7, -1.5, 5.5);

  FailingProbePt        = fs->make<TH1F>("FailingProbePt", "FailingProbe pT", 10000, 0, 10000);
  FailingProbeEta       = fs->make<TH1F>("FailingProbeEta", "FailingProbe #eta",    96, -4.8, 4.8);
  FailingProbePhi       = fs->make<TH1F>("FailingProbePhi", "FailingProbe #phi", 41, -TMath::Pi(), TMath::Pi());
  FailingProbeNVertices = fs->make<TH1F>("FailingProbeNVertices", "FailingProbe # vertices/event",  200, 0, 200);
  FailingProbeEtaPhi    = fs->make<TH2F>("FailingProbeEtaPhi", "FailingProbe #eta #phi",    14, eta_bins_for_2D, 41, -TMath::Pi(), TMath::Pi());
  FailingProbeEtaPt     = fs->make<TH2F>("FailingProbeEtaPt", "FailingProbe #eta pT",    14, eta_bins_for_2D, 10000, 0, 10000);
  FailingProbeEtaDptPt  = fs->make<TH2F>("FailingProbeEtaDptPt", "FailingProbe #eta dpt/pT",    14, eta_bins_for_2D, 1000, 0, 2);
  if(isAOD)
    FailingProbeEtaShower = fs->make<TH2F>("FailingProbeEtaShower", "FailingProbe #eta shower",    14, eta_bins_for_2D, 7, -1.5, 5.5);

  // TnP pair
  PairNoPtMass    = fs->make<TH1F>("PairNoPtMass", "TnP PairNoPt mass", 10000, 0, 10000);
  PairNoPtPt      = fs->make<TH1F>("PairNoPtPt", "TnP PairNoPt pT", 10000, 0, 10000);
  PairNoPtEta     = fs->make<TH1F>("PairNoPtEta", "TnP PairNoPt #eta",    96, -4.8, 4.8);
  PairNoPtRap     = fs->make<TH1F>("PairNoPtRap", "TnP PairNoPt y", 96, -4.8, 4.8);

  PairMass        = fs->make<TH1F>("PairMass", "TnP Pair mass", 10000, 0, 10000);
  PairPt          = fs->make<TH1F>("PairPt", "TnP Pair pT", 10000, 0, 10000);
  PairEta         = fs->make<TH1F>("PairEta", "TnP Pair #eta",    96, -4.8, 4.8);
  PairRap         = fs->make<TH1F>("PairRap", "TnP Pair y", 96, -4.8, 4.8);

  PassingPairMass = fs->make<TH1F>("PassingPairMass", "TnP PassingPair mass", 10000, 0, 10000);
  PassingPairPt   = fs->make<TH1F>("PassingPairPt", "TnP PassingPair pT", 10000, 0, 10000);
  PassingPairEta  = fs->make<TH1F>("PassingPairEta", "TnP PassingPair #eta",    96, -4.8, 4.8);
  PassingPairRap  = fs->make<TH1F>("PassingPairRap", "TnP PassingPair y", 96, -4.8, 4.8);

  FailingPairMass = fs->make<TH1F>("FailingPairMass", "TnP FailingPair mass", 10000, 0, 10000);
  FailingPairPt   = fs->make<TH1F>("FailingPairPt", "TnP FailingPair pT", 10000, 0, 10000);
  FailingPairEta  = fs->make<TH1F>("FailingPairEta", "TnP FailingPair #eta",    96, -4.8, 4.8);
  FailingPairRap  = fs->make<TH1F>("FailingPairRap", "TnP FailingPair y", 96, -4.8, 4.8);

  // Probe ID variable
  /*
    ProbeAbsTkIso = fs->make<TH1F>("ProbeAbsTkIso", "Probe Iso. (#Delta R < 0.3) #Sigma pT", 1000, 0, 1000);
    ProbeRelTkIso = fs->make<TH1F>("ProbeRelTkIso", "Probe Iso. (#Delta R < 0.3) #Sigma pT / tk. pT", 500, 0, 5);
    ProbeChi2dof = fs->make<TH1F>("ProbeChi2dof", "Probe #chi^{2}/dof", 500, 0, 50);
    ProbeTrackDXYBS = fs->make<TH1F>("ProbeTrackDXYBS", "Probe |dxy wrt BS|", 10000, 0, 2);
    ProbeTrackDZBS = fs->make<TH1F>("ProbeTrackDZBS", "Probe |dz wrt BS|", 10000, 0, 20);
    ProbeTrackDXYPV = fs->make<TH1F>("ProbeTrackDXYPV", "Probe |dxy wrt PV|", 10000, 0, 2);
    ProbeTrackDZPV = fs->make<TH1F>("ProbeTrackDZPV", "Probe |dz wrt PV|", 10000, 0, 20);
    ProbeNPxHits = fs->make<TH1F>("ProbeNPxHits", "Probe # pixel hits", 10, 0,  10);
    ProbeNStHits = fs->make<TH1F>("ProbeNStHits", "Probe # strip hits", 40, 0, 40);
    ProbeNTkHits = fs->make<TH1F>("ProbeNTkHits", "Probe # tracker hits", 50, 0, 50);
    ProbeNMuHits = fs->make<TH1F>("ProbeNMuHits", "Probe # muon hits", 60, 0, 60);
    ProbeNHits = fs->make<TH1F>("ProbeNHits", "Probe # hits", 80, 0, 80);
    ProbeNPxLayers = fs->make<TH1F>("ProbeNPxLayers", "Probe # pixel layers", 10, 0, 10);
    ProbeNStLayers = fs->make<TH1F>("ProbeNStLayers", "Probe # strip layers", 20, 0, 20);
    ProbeNTkLayers = fs->make<TH1F>("ProbeNTkLayers", "Probe # tracker layers", 30, 0, 30);
  */

  // Bin histograms
  if(useBinHistos2D) {
    Pt2D      = make_bin_histos_2D("Pt",  vec_PtBins);
    PtB2D     = make_bin_histos_2D("PtB", vec_PtBins);
    PtO2D     = make_bin_histos_2D("PtO", vec_PtBins);
    PtE2D     = make_bin_histos_2D("PtE", vec_PtBins);
    PtF2D     = make_bin_histos_2D("PtF", vec_PtBins);

    AbsP2D    = make_bin_histos_2D("AbsP",  vec_AbsPBins);
    AbsPB2D   = make_bin_histos_2D("AbsPB", vec_AbsPBins);
    AbsPO2D   = make_bin_histos_2D("AbsPO", vec_AbsPBins);
    AbsPE2D   = make_bin_histos_2D("AbsPE", vec_AbsPBins);
    AbsPF2D   = make_bin_histos_2D("AbsPF", vec_AbsPBins);

    Eta2D     = make_bin_histos_2D("Eta", vec_EtaBins);
    Phi2D     = make_bin_histos_2D("Phi", vec_PhiBins);
    Vtx2D     = make_bin_histos_2D("Vtx", vec_VtxBins);

    if(isAOD) {
      Shower2D  = make_bin_histos_2D("Shower",  vec_ShowerBins);
      ShowerB2D = make_bin_histos_2D("ShowerB", vec_ShowerBins);
      ShowerE2D = make_bin_histos_2D("ShowerE", vec_ShowerBins);
    }
  }

  else {
    Pt      = make_bin_histos("Pt",  vec_PtBins);
    PtB     = make_bin_histos("PtB", vec_PtBins);
    PtO     = make_bin_histos("PtO", vec_PtBins);
    PtE     = make_bin_histos("PtE", vec_PtBins);
    PtF     = make_bin_histos("PtF", vec_PtBins);

    AbsP    = make_bin_histos("AbsP",  vec_AbsPBins);
    AbsPB   = make_bin_histos("AbsPB", vec_AbsPBins);
    AbsPO   = make_bin_histos("AbsPO", vec_AbsPBins);
    AbsPE   = make_bin_histos("AbsPE", vec_AbsPBins);
    AbsPF   = make_bin_histos("AbsPF", vec_AbsPBins);

    Eta     = make_bin_histos("Eta", vec_EtaBins);
    Phi     = make_bin_histos("Phi", vec_PhiBins);
    Vtx     = make_bin_histos("Vtx", vec_VtxBins);
    if(isAOD) {
      Shower  = make_bin_histos("Shower",  vec_ShowerBins);
      ShowerB = make_bin_histos("ShowerB", vec_ShowerBins);
      ShowerE = make_bin_histos("ShowerE", vec_ShowerBins);
    }
  }

  // Templet histograms
  double *arr_PtBins;
  double *arr_AbsPBins;
  double *arr_EtaBins;
  double *arr_PhiBins;
  double *arr_VtxBins;
  double *arr_ShowerBins;
  CopyVectorToArray(vec_PtBins, arr_PtBins);
  CopyVectorToArray(vec_AbsPBins, arr_AbsPBins);
  CopyVectorToArray(vec_EtaBins, arr_EtaBins);
  CopyVectorToArray(vec_PhiBins, arr_PhiBins);
  CopyVectorToArray(vec_VtxBins, arr_VtxBins);
  CopyVectorToArray(vec_ShowerBins, arr_ShowerBins);
  hEffTemplatePt     = fs->make<TH1F>("hEffTemplatePt", "", (int)vec_PtBins.size()-1, arr_PtBins);
  hEffTemplateAbsP   = fs->make<TH1F>("hEffTemplateAbsP", "", (int)vec_AbsPBins.size()-1, arr_AbsPBins);
  hEffTemplateEta    = fs->make<TH1F>("hEffTemplateEta", "", (int)vec_EtaBins.size()-1, arr_EtaBins);
  hEffTemplatePhi    = fs->make<TH1F>("hEffTemplatePhi", "", (int)vec_PhiBins.size()-1, arr_PhiBins);
  hEffTemplateVtx    = fs->make<TH1F>("hEffTemplateVtx", "", (int)vec_VtxBins.size()-1, arr_VtxBins);
  hEffTemplateShower = fs->make<TH1F>("hEffTemplateShower", "", (int)vec_ShowerBins.size()-1, arr_ShowerBins);

  comparison_tree = fs->make<TTree>("t", "");
  comparison_tree->Branch("IsRealData", &IsRealData, "IsRealData/O");
  comparison_tree->Branch("RunNum",&RunNum,"RunNum/I");
  comparison_tree->Branch("LumiBlockNum",&LumiBlockNum,"LumiBlockNum/I");
  comparison_tree->Branch("EventNum",&EventNum,"EventNum/l");
  comparison_tree->Branch("Mass",&Mass,"Mass/D");
  comparison_tree->Branch("VertexMass",&VertexMass,"VertexMass/D");
  comparison_tree->Branch("Probe_Pt",&Probe_Pt,"Probe_Pt/D");
  comparison_tree->Branch("Probe_Eta",&Probe_Eta,"Probe_Eta/D");
}

void CustoTnPHistosForTnP::getBSandPV(const edm::Event& event) {
  // We store these as bare pointers. Should find better way, but
  // don't want to pass them around everywhere...
  edm::Handle<reco::BeamSpot> hbs;
  event.getByLabel(beamspot_src, hbs);
  beamspot = hbs.isValid() ? &*hbs : 0; // nice and fragile
  NBeamSpot->Fill(beamspot != 0);

  edm::Handle<reco::VertexCollection> vertices;
  event.getByLabel(vertex_src, vertices);
  vertex = 0;
  int vertex_count = 0;
  for (reco::VertexCollection::const_iterator it = vertices->begin(), ite = vertices->end(); it != ite; ++it) {
    if (it->ndof() > 4 && fabs(it->z()) <= 24 && fabs(it->position().rho()) <= 2) {
      if (vertex == 0)
        vertex = &*it;
      ++vertex_count;
    }
  }
  nVtx = vertex_count;
  NVertices->Fill(vertex_count, _totalWeight );
}

int CustoTnPHistosForTnP::calcNShowers(const reco::CandidateBaseRef& mu) {

  int etaCat = -1;
  if( fabs(mu->eta()) < 0.9 )
    etaCat = 1;
  else if( fabs(mu->eta()) >= 1.2 )
    etaCat = 2;

  if( etaCat<0 )
    return -999;

  int nShowers = 0;

  const pat::Muon* muPat = toConcretePtr<pat::Muon>(mu);
  int st1 = muPat->userInt("nHits1");
  int st2 = muPat->userInt("nHits2");
  int st3 = muPat->userInt("nHits3");
  int st4 = muPat->userInt("nHits4");

  if(etaCat==1) {
    if(st1>26)
      nShowers +=1;
    if(st2>26)
      nShowers +=1;
    if(st3>26)
      nShowers +=1;
    // if(st4>18)
    //   nShowers +=1;
  }

  else if(etaCat==2) {
    if(st1>18)
      nShowers +=1;
    if(st2>18)
      nShowers +=1;
    if(st3>18)
      nShowers +=1;
    if(st4>18)
      nShowers +=1;
  }

  std::cout << "\netaCat: " << etaCat << std::endl;
  std::cout << "\tst1= " << st1 << std::endl;
  std::cout << "\tst2= " << st2 << std::endl;
  std::cout << "\tst3= " << st3 << std::endl;
  std::cout << "\tst4= " << st4 << std::endl;
  std::cout << "\t--> nShowers= " << nShowers << std::endl;

  return nShowers;
}

void CustoTnPHistosForTnP::fillTnPControlHistos(const pat::CompositeCandidate& dil,
                                                 const reco::CandidateBaseRef& TagMu,
                                                 const reco::CandidateBaseRef& ProbeMu,
                                                 float probe_dpt_over_pt,
                                                 int   probe_nshowers,
                                                 bool  isPassNoPt ) {

  //for offline variables
  // const pat::Muon* TagPat = toConcretePtr<pat::Muon>(TagMu);
  // const pat::Muon* ProbePat = toConcretePtr<pat::Muon>(ProbeMu);

  TagPt->Fill( TagMu->pt(), _totalWeight );
  TagEta->Fill( TagMu->eta(), _totalWeight );
  TagPhi->Fill( TagMu->phi(), _totalWeight );

  /*
  if(TagPat) {  // fill offline
    const reco::MuonIsolation& iso = TagPat->isolationR03();
    TagAbsTkIso->Fill( iso.sumPt, _totalWeight );
    TagRelTkIso->Fill( iso.sumPt / TagPat->innerTrack()->pt(), _totalWeight );

    //const reco::TrackRef track = patmuon::getPickedTrack(*TagPat);
    const reco::TrackRef GlTrack = TagPat->globalTrack();
    const reco::TrackRef InTrack = TagPat->innerTrack();
    const reco::TrackRef BeTrack = TagPat->muonBestTrack();
    if (GlTrack.isAvailable() && InTrack.isAvailable() && BeTrack.isAvailable()) {
      TagChi2dof->Fill( GlTrack->normalizedChi2(), _totalWeight );

      if (beamspot != 0) {
        TagTrackDXYBS->Fill( fabs(BeTrack->dxy(beamspot->position())), _totalWeight );
        TagTrackDZBS->Fill( fabs(BeTrack->dz (beamspot->position())), _totalWeight );
      }

      if (vertex != 0) {
        TagTrackDXYPV->Fill( fabs(BeTrack->dxy(vertex->position())), _totalWeight );
        TagTrackDZPV->Fill( fabs(BeTrack->dz (vertex->position())), _totalWeight );
      }

      const reco::HitPattern& hp = InTrack->hitPattern();
      TagNPxHits->Fill( hp.numberOfValidPixelHits(), _totalWeight );
      TagNStHits->Fill( hp.numberOfValidStripHits(), _totalWeight );
      TagNTkHits->Fill( hp.numberOfValidTrackerHits(), _totalWeight );
      TagNMuHits->Fill( GlTrack->hitPattern().numberOfValidMuonHits(), _totalWeight );

      TagNHits->Fill( GlTrack->hitPattern().numberOfValidHits(), _totalWeight );

      TagNPxLayers->Fill( hp.pixelLayersWithMeasurement(), _totalWeight );
      TagNStLayers->Fill( hp.stripLayersWithMeasurement(), _totalWeight );
      TagNTkLayers->Fill( hp.trackerLayersWithMeasurement(), _totalWeight );
    }
  }
  */

  ProbePt->Fill( ProbeMu->pt(), _totalWeight );
  if( isPassNoPt )
    PassingProbePt->Fill( ProbeMu->pt(), _totalWeight );
  else
    FailingProbePt->Fill( ProbeMu->pt(), _totalWeight );

  PairNoPtMass->Fill( dil.mass(), _totalWeight );
  PairNoPtPt->Fill( dil.pt(), _totalWeight );
  PairNoPtEta->Fill( dil.eta(), _totalWeight );
  PairNoPtRap->Fill( dil.rapidity(), _totalWeight );

  if(ProbeMu->pt() > probe_pt_min) {
    ProbeEta->Fill( ProbeMu->eta(), _totalWeight );
    ProbePhi->Fill( ProbeMu->phi(), _totalWeight );
    ProbeNVertices->Fill( nVtx, _totalWeight );
    ProbeEtaPhi->Fill( ProbeMu->eta(), ProbeMu->phi(), _totalWeight );
    ProbeEtaPt->Fill( ProbeMu->eta(), ProbeMu->pt(), _totalWeight );
    ProbeEtaDptPt->Fill( ProbeMu->eta(), probe_dpt_over_pt, _totalWeight );
    if(isAOD)
      ProbeEtaShower->Fill( ProbeMu->eta(), probe_nshowers, _totalWeight );

    PairMass->Fill( dil.mass(), _totalWeight );
    PairPt->Fill( dil.pt(), _totalWeight );
    PairEta->Fill( dil.eta(), _totalWeight );
    PairRap->Fill( dil.rapidity(), _totalWeight );

    if( isPassNoPt ) {
      PassingProbeEta->Fill( ProbeMu->eta(), _totalWeight );
      PassingProbePhi->Fill( ProbeMu->phi(), _totalWeight );
      PassingProbeNVertices->Fill( nVtx, _totalWeight );
      PassingProbeEtaPhi->Fill( ProbeMu->eta(), ProbeMu->phi(), _totalWeight );
      PassingProbeEtaPt->Fill( ProbeMu->eta(), ProbeMu->pt(), _totalWeight );
      PassingProbeEtaDptPt->Fill( ProbeMu->eta(), probe_dpt_over_pt, _totalWeight );
      if(isAOD)
        PassingProbeEtaShower->Fill( ProbeMu->eta(), probe_nshowers, _totalWeight );

      PassingPairMass->Fill( dil.mass(), _totalWeight );
      PassingPairPt->Fill( dil.pt(), _totalWeight );
      PassingPairEta->Fill( dil.eta(), _totalWeight );
      PassingPairRap->Fill( dil.rapidity(), _totalWeight );
    }
    else {
      FailingProbeEta->Fill( ProbeMu->eta(), _totalWeight );
      FailingProbePhi->Fill( ProbeMu->phi(), _totalWeight );
      FailingProbeNVertices->Fill( nVtx, _totalWeight );
      FailingProbeEtaPhi->Fill( ProbeMu->eta(), ProbeMu->phi(), _totalWeight );
      FailingProbeEtaPt->Fill( ProbeMu->eta(), ProbeMu->pt(), _totalWeight );
      FailingProbeEtaDptPt->Fill( ProbeMu->eta(), probe_dpt_over_pt, _totalWeight );
      if(isAOD)
        FailingProbeEtaShower->Fill( ProbeMu->eta(), probe_nshowers, _totalWeight );

      FailingPairMass->Fill( dil.mass(), _totalWeight );
      FailingPairPt->Fill( dil.pt(), _totalWeight );
      FailingPairEta->Fill( dil.eta(), _totalWeight );
      FailingPairRap->Fill( dil.rapidity(), _totalWeight );
    }

    /*
    if(ProbePat) {  // fill offline
      const reco::MuonIsolation& iso = ProbePat->isolationR03();
      ProbeAbsTkIso->Fill( iso.sumPt, _totalWeight );
      ProbeRelTkIso->Fill( iso.sumPt / ProbePat->innerTrack()->pt(), _totalWeight );

      //const reco::TrackRef track = patmuon::getPickedTrack(*ProbePat);
      const reco::TrackRef GlTrack = ProbePat->globalTrack();
      const reco::TrackRef InTrack = ProbePat->innerTrack();
      const reco::TrackRef BeTrack = ProbePat->muonBestTrack();
      if (GlTrack.isAvailable() && InTrack.isAvailable() && BeTrack.isAvailable()) {
        ProbeChi2dof->Fill( GlTrack->normalizedChi2(), _totalWeight );

        if (beamspot != 0) {
          ProbeTrackDXYBS->Fill( fabs(BeTrack->dxy(beamspot->position())), _totalWeight );
          ProbeTrackDZBS->Fill( fabs(BeTrack->dz (beamspot->position())), _totalWeight );
        }

        if (vertex != 0) {
          ProbeTrackDXYPV->Fill( fabs(BeTrack->dxy(vertex->position())), _totalWeight );
          ProbeTrackDZPV->Fill( fabs(BeTrack->dz (vertex->position())), _totalWeight );
        }

        const reco::HitPattern& hp = InTrack->hitPattern();
        ProbeNPxHits->Fill( hp.numberOfValidPixelHits(), _totalWeight );
        ProbeNStHits->Fill( hp.numberOfValidStripHits(), _totalWeight );
        ProbeNTkHits->Fill( hp.numberOfValidTrackerHits(), _totalWeight );
        ProbeNMuHits->Fill( GlTrack->hitPattern().numberOfValidMuonHits(), _totalWeight );

        ProbeNHits->Fill( GlTrack->hitPattern().numberOfValidHits(), _totalWeight );

        ProbeNPxLayers->Fill( hp.pixelLayersWithMeasurement(), _totalWeight );
        ProbeNStLayers->Fill( hp.stripLayersWithMeasurement(), _totalWeight );
        ProbeNTkLayers->Fill( hp.trackerLayersWithMeasurement(), _totalWeight );
      }
    }
    */
  }
}

void CustoTnPHistosForTnP::fillTnPBinHistos(  double               dil_mass,
                                              double               probe_pt,
                                              bool                 isPassNoPt,
                                              double               binValue,
                                              std::vector<double>& vec_bins,
                                              BinHistos&           Histos,
                                              bool hasPtCut = true ) {

  if( hasPtCut && (probe_pt < probe_pt_min) )
    return;

  bool isPass = false;

  if( hasPtCut )
    isPass = ( isPassNoPt && (probe_pt > probe_pt_min) );
  else
    isPass = isPassNoPt;

  int nbins = ( (int)vec_bins.size() ) - 1;

  for(int i=0; i<nbins; ++i) {
    if( binValue >= vec_bins[i] && binValue < vec_bins[i+1] ) {

      if(isPass)
        (Histos.first)[i]->Fill( dil_mass, _totalWeight );
      else
        (Histos.second)[i]->Fill( dil_mass, _totalWeight );

      break;
    }
  }

}

void CustoTnPHistosForTnP::fillTnPBinHistos2D(  double               dil_mass,
                                                double               probe_pt,
                                                bool                 isPassNoPt,
                                                double               binValue,
                                                BinHistos2D&         Histos,
                                                bool hasPtCut = true ) {

  if( hasPtCut && (probe_pt < probe_pt_min) )
    return;

  bool isPass = false;

  if( hasPtCut )
    isPass = ( isPassNoPt && (probe_pt > probe_pt_min) );
  else
    isPass = isPassNoPt;

  if(isPass)
    Histos[0]->Fill( dil_mass, binValue, _totalWeight );
  else
    Histos[1]->Fill( dil_mass, binValue, _totalWeight );

  Histos[2]->Fill( dil_mass, binValue, binValue*_totalWeight );
  Histos[3]->Fill( dil_mass, binValue, binValue*binValue*_totalWeight );

}


void CustoTnPHistosForTnP::analyze(const edm::Event& event, const edm::EventSetup& setup) {

  IsRealData = false;
  RunNum = -999;
  LumiBlockNum = -999;
  EventNum = 0;
  Mass = -999.;
  VertexMass = -999.;
  Probe_Pt = -999.;
  Probe_Eta = -999.;

  //---- Prescales : Not using for now...
    //  edm::Handle<int> hltPrescale;
    //  edm::Handle<int> l1Prescale;
    //  event.getByLabel(edm::InputTag("getPrescales","HLTPrescale","CustoTnPAnalysis"), hltPrescale);
    //  event.getByLabel(edm::InputTag("getPrescales","L1Prescale","CustoTnPAnalysis"), l1Prescale);
    if (_usePrescaleWeight) {
      edm::Handle<int> totalPrescale;
      event.getByLabel(edm::InputTag("getPrescales","TotalPrescale","CustoTnPAnalysis"), totalPrescale);
      _prescaleWeight = *totalPrescale;
    }
    //  std::cout<<*hltPrescale<<std::endl;
    //  std::cout<<l1Prescale<<std::endl;
    //  std::cout<<totalPrescale<<std::endl;

  //---- Generator weights
  if (_useMadgraphWeight) {
    _eventWeight = 1.;
    _madgraphWeight = 1.;

    edm::Handle<GenEventInfoProduct> gen_ev_info;
    event.getByLabel(edm::InputTag("generator"), gen_ev_info);
    if (gen_ev_info.isValid()){
      _eventWeight = gen_ev_info->weight();
      _madgraphWeight = ( _eventWeight > 0 ) ? 1.0 : -1.0;
    }
    WeightMadGraph->Fill( _madgraphWeight );
  }

  //---- Get PileUp Weights
  int thePU = -1;
  _pileupWeight = 1.;
  if( !event.isRealData() ){

    edm::Handle<std::vector< PileupSummaryInfo > > pileup;
    if ( event.getByLabel(pileup_src, pileup)){
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = pileup->begin(); PVI != pileup->end(); ++PVI)
      {
        if(PVI->getBunchCrossing()==0){
          thePU = PVI->getTrueNumInteractions();
          continue;
        }
      }

      int nPileUpBin = (int)vec_PileUpWeight.size();
      for(int iPU=0; iPU<nPileUpBin; ++iPU) {
        if(thePU == iPU) {
          _pileupWeight = vec_PileUpWeight[iPU];
        }
      }
    }
    else
      edm::LogError("") << "PU collection not found !!!";

    if(!ShutUp)  std::cout << "CustoTnPHistosForTnP::analyze : PU = " << thePU << "  PU weight = " << _pileupWeight << std::endl;
  }

  _totalWeight = 1.;
  if( _madgraphWeight != 1. || _pileupWeight != 1. )
    _totalWeight = _madgraphWeight * _pileupWeight;  // prescale weight is not considered
  if(!ShutUp)  std::cout << "CustoTnPHistosForTnP::analyze : Total weight = " << _totalWeight << std::endl;

  if (use_bs_and_pv)
    getBSandPV(event);

  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  event.getByLabel(dilepton_src, dileptons);

  if( !dileptons.isValid() ) {
    std::cout << "CustoTnPHistosForTnP::analyze : !dileptons.isValid() ---> return" << std::endl;
    return;
  }

  NDileptons->Fill(dileptons->size(), _totalWeight );

  pat::CompositeCandidateCollection::const_iterator dil = dileptons->begin(), dile = dileptons->end();
  for ( ; dil != dile; ++dil) {

    if( !( dil->hasUserInt("isTag0Probe1") && dil->hasUserInt("isTag1Probe0") ) ) {
      std::cout << "CustoTnPHistosForTnP::analyze : no isTag0Probe1 or isTag1Probe0 in dil" << std::endl;
      continue;
    }

    double dil_mass = dil->mass();
    double dil_vertex_mass = dil->userFloat("vertexM");
    float lep0_dpt_over_pt = dil->userFloat("lep0_dpt_over_pt");
    float lep1_dpt_over_pt = dil->userFloat("lep1_dpt_over_pt");

    const reco::CandidateBaseRef& lep0 = dileptonDaughter(*dil, 0);
    const reco::CandidateBaseRef& lep1 = dileptonDaughter(*dil, 1);

    if(lep0.isNonnull() && lep1.isNonnull()) {
      const pat::Muon* mu0 = toConcretePtr<pat::Muon>(lep0);
      const pat::Muon* mu1 = toConcretePtr<pat::Muon>(lep1);
      if(mu0 && mu1) {

        bool isAleadyFilled = false;

        //---- Tag0 and Probe1
        if( dil->userInt("isTag0Probe1") ) {
          const reco::CandidateBaseRef& TagMu = lep0;
          const reco::CandidateBaseRef& ProbeMu = lep1;

          int probe_nshowers = -999;
          if(isAOD)
            probe_nshowers = calcNShowers(ProbeMu);

          if(!ShutUp)  std::cout << "CustoTnPHistosForTnP::analyze : Tag0 and Probe1" << std::endl;
          if(!ShutUp)  std::cout << "                                              pT=" << ProbeMu->pt() << std::endl;
          if(!ShutUp)  std::cout << "                                             eta=" << ProbeMu->eta() << std::endl;
          if(!ShutUp)  std::cout << "                                             phi=" << ProbeMu->phi() << std::endl;
          if(!ShutUp)  std::cout << "                                            nVtx=" << nVtx << std::endl;

          bool isPassingProbe = ( passing_probe_selector(*mu1) && (lep1_dpt_over_pt < passing_probe_dpt_over_pt_max) && (fabs(mu1->innerTrack()->dz( vertex->position() )) < passing_probe_dz_max) );
          if(!ShutUp)  std::cout << "                                  isPassingProbe = " << isPassingProbe << std::endl;

          fillTnPControlHistos(*dil, TagMu, ProbeMu, lep1_dpt_over_pt, probe_nshowers, isPassingProbe);

          if(useBinHistos2D) {
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), Pt2D, false);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtB2D, false);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtO2D, false);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtE2D, false);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtF2D, false);

            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsP2D, true);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPB2D, true);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPO2D, true);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPE2D, true);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPF2D, true);

            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->eta(), Eta2D, true);
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->phi(), Phi2D, true);
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, nVtx, Vtx2D, true);

            if(isAOD) {
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, Shower2D, true);
              if( fabs(ProbeMu->eta())<0.9 )
                fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerB2D, true);
              else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
                fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerE2D, true);
            }
          }
          else {
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, Pt, false);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtB, false);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtO, false);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtE, false);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtF, false);

            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsP, true);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPB, true);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPO, true);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPE, true);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPF, true);

            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->eta(), vec_EtaBins, Eta, true);
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->phi(), vec_PhiBins, Phi, true);
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, nVtx, vec_VtxBins, Vtx, true);

            if(isAOD) {
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, Shower, true);
              if( fabs(ProbeMu->eta())<0.9 )
                fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerB, true);
              else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
                fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerE, true);
            }
          }

          if( ProbeMu->pt() > probe_pt_min ) {
            bool isComparisonProbe = ( comparison_probe_selector(*mu1) && (lep1_dpt_over_pt < comparison_probe_dpt_over_pt_max) && (fabs(mu1->innerTrack()->dz( vertex->position() )) < comparison_probe_dz_max) );
            if(!ShutUp)  std::cout << "                               isComparisonProbe = " << isComparisonProbe << std::endl;
            if( isPassingProbe && !isComparisonProbe ) {
              IsRealData = event.isRealData();
              RunNum = event.id().run();
              LumiBlockNum = event.id().luminosityBlock();
              EventNum = event.id().event();
              Mass = dil_mass;
              VertexMass = dil_vertex_mass;
              Probe_Pt = ProbeMu->pt();
              Probe_Eta = ProbeMu->eta();
              comparison_tree->Fill();
              isAleadyFilled = true;
            }
          }
        } // Tag0 and Probe1

        //---- Tag1 and Probe0
        if( dil->userInt("isTag1Probe0") ) {
          const reco::CandidateBaseRef& TagMu = lep1;
          const reco::CandidateBaseRef& ProbeMu = lep0;

          int probe_nshowers = -999;
          if(isAOD)
            probe_nshowers = calcNShowers(ProbeMu);

          if(!ShutUp)  std::cout << "CustoTnPHistosForTnP::analyze : Tag1 and Probe0" << std::endl;
          if(!ShutUp)  std::cout << "                                              pT=" << ProbeMu->pt() << std::endl;
          if(!ShutUp)  std::cout << "                                             eta=" << ProbeMu->eta() << std::endl;
          if(!ShutUp)  std::cout << "                                             phi=" << ProbeMu->phi() << std::endl;
          if(!ShutUp)  std::cout << "                                            nVtx=" << nVtx << std::endl;

          bool isPassingProbe = ( passing_probe_selector(*mu0) && (lep0_dpt_over_pt < passing_probe_dpt_over_pt_max) && (fabs(mu0->innerTrack()->dz( vertex->position() )) < passing_probe_dz_max) );
          if(!ShutUp)  std::cout << "                                  isPassingProbe = " << isPassingProbe << std::endl;

          fillTnPControlHistos(*dil, TagMu, ProbeMu, lep0_dpt_over_pt, probe_nshowers, isPassingProbe);

          if(useBinHistos2D) {
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), Pt2D, false);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtB2D, false);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtO2D, false);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtE2D, false);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), PtF2D, false);

            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsP2D, true);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPB2D, true);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPO2D, true);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPE2D, true);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), AbsPF2D, true);

            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->eta(), Eta2D, true);
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->phi(), Phi2D, true);
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, nVtx, Vtx2D, true);

            if(isAOD) {
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, Shower2D, true);
              if( fabs(ProbeMu->eta())<0.9 )
                fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerB2D, true);
              else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
                fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, ShowerE2D, true);
            }
          }
          else {
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, Pt, false);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtB, false);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtO, false);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtE, false);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->pt(), vec_PtBins, PtF, false);

            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsP, true);
            if( fabs(ProbeMu->eta())<0.9 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPB, true);
            else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPO, true);
            else if( fabs(ProbeMu->eta())>=1.2 && fabs(ProbeMu->eta())<2.1 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPE, true);
            else if( fabs(ProbeMu->eta())>=2.1 && fabs(ProbeMu->eta())<2.4 )
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->p(), vec_AbsPBins, AbsPF, true);

            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->eta(), vec_EtaBins, Eta, true);
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, ProbeMu->phi(), vec_PhiBins, Phi, true);
            fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, nVtx, vec_VtxBins, Vtx, true);

            if(isAOD) {
              fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, Shower, true);
              if( fabs(ProbeMu->eta())<0.9 )
                fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerB, true);
              else if( fabs(ProbeMu->eta())>=0.9 && fabs(ProbeMu->eta())<1.2 )
                fillTnPBinHistos(dil_mass, ProbeMu->pt(), isPassingProbe, probe_nshowers, vec_ShowerBins, ShowerE, true);
            }
          }

          if( ProbeMu->pt() > probe_pt_min ) {
            bool isComparisonProbe = ( comparison_probe_selector(*mu0) && (lep0_dpt_over_pt < comparison_probe_dpt_over_pt_max) && (fabs(mu0->innerTrack()->dz( vertex->position() )) < comparison_probe_dz_max) );
            if(!ShutUp)  std::cout << "                               isComparisonProbe = " << isComparisonProbe << std::endl;
            if( !isAleadyFilled && isPassingProbe && !isComparisonProbe ) {
              IsRealData = event.isRealData();
              RunNum = event.id().run();
              LumiBlockNum = event.id().luminosityBlock();
              EventNum = event.id().event();
              Mass = dil_mass;
              VertexMass = dil_vertex_mass;
              Probe_Pt = ProbeMu->pt();
              Probe_Eta = ProbeMu->eta();
              comparison_tree->Fill();
            }
          }
        } // Tag1 and Probe0

        if(!ShutUp)  std::cout << std::endl;
      } // if(mu0 && mu1)
    } // if(lep0.isNonnull() && lep1.isNonnull())


  } // for ( ; dil != dile; ++dil)

}

DEFINE_FWK_MODULE(CustoTnPHistosForTnP);
