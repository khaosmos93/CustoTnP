// -*- C++ -*-
//
// Package:    Histograms
// Class:      EventCounter
// 
/**\class EventCounter EventCounter.cc brot/EventCounter/src/EventCounter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>

**/

// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <iostream> 

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


//ROOT
#include "TH1.h"
//
// class decleration
//

class EventCounter : public edm::EDAnalyzer {
public:
  explicit EventCounter(const edm::ParameterSet&);
  ~EventCounter();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::EDGetTokenT<GenEventInfoProduct> GenInfoTag_;

  TH1F* Events;
  TH1F* weights;
};

// constructors and destructor
EventCounter::EventCounter(const edm::ParameterSet& iConfig):
  GenInfoTag_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfoTag")))
{
  edm::Service<TFileService> fs;
  TH1::SetDefaultSumw2(true);

  Events  = fs->make<TH1F>("Events",  "Events",  4, -2, 2);
  weights = fs->make<TH1F>("weights", "weights", 4, -2, 2);
}



EventCounter::~EventCounter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// member functions
// ------------ method called to for each event  ------------
void
EventCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  Events->Fill(1.0);

  //get generator weigts
  if( !iEvent.isRealData() ) {
    edm::Handle<GenEventInfoProduct> genInfoProduct;
    iEvent.getByToken(GenInfoTag_, genInfoProduct);

    if( genInfoProduct.isValid() ) {
      if ((*genInfoProduct).weight() < 0.0) {
        weights->Fill(-1.0);
      }
      else {
        weights->Fill(1.0);
      }
    }
    else {
      weights->Fill(0.0);
    }
  }
  else {
    weights->Fill(1.0);
  }

}

// ------------ method called once each job just before starting event loop  ------------
void 
EventCounter::beginJob()
{
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EventCounter::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventCounter);
