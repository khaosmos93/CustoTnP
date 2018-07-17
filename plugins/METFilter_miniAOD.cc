// -*- C++ -*-
//
// Package:    METFilter_miniAOD
// Class:      METFilter_miniAOD
// 
/**\class METFilter_miniAOD METFilter_miniAOD.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/



// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"


#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
//
// class declaration
//

class METFilter_miniAOD : public edm::EDFilter {
public:
  explicit METFilter_miniAOD(const edm::ParameterSet&);
  ~METFilter_miniAOD();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  
  edm::EDGetTokenT<edm::TriggerResults> filterToken_;
  std::string filterFlag_;



  bool debug;
};

// constructors and destructor
METFilter_miniAOD::METFilter_miniAOD(const edm::ParameterSet& iConfig):
  filterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("src")))
{
  filterFlag_    = iConfig.getParameter < std::string > ("flag");

  debug = false;
}

METFilter_miniAOD::~METFilter_miniAOD(){}


// member functions
// ------------ method called on each new Event  ------------
bool
METFilter_miniAOD::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::TriggerResults> filterBits;
  iEvent.getByToken(filterToken_, filterBits);
  const edm::TriggerNames &allFilterNames = iEvent.triggerNames(*filterBits);
  bool filter = false;
  
  for (unsigned int i = 0, n = filterBits->size(); i < n; ++i){
	  if (allFilterNames.triggerName(i) == filterFlag_ && filterBits->accept(i)){
		  filter = true;
	  }
  }

  return filter;
}

// ------------ method called once each job just before starting event loop  ------------
void 
METFilter_miniAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METFilter_miniAOD::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(METFilter_miniAOD);
