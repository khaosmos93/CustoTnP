import FWCore.ParameterSet.Config as cms

# A filter for post-tuple filtering on the goodData results as stored
# in a TriggerResults object instead of filtering at tuple-making
# time.
import HLTrigger.HLTfilters.hltHighLevel_cfi
goodDataFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
goodDataFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', 'PAT')
goodDataFilter.HLTPaths = ["goodDataPrimaryVertexFilter"] # can set to just 'goodDataPrimaryVertexFilter', for example
#goodDataFilter.HLTPaths = ['goodDataMETFilter']
goodDataFilter.andOr = False # = AND

from CustoTnP.Analyzer.goodData_cff import primaryVertexMiniAOD, hltPhysicsDeclared, metFilters
goodDataFiltersMiniAOD = [primaryVertexMiniAOD]
## for full filtering, use:
#goodDataFiltersMiniAOD = [primaryVertexMiniAOD,hltPhysicsDeclared]
#goodDataFiltersMiniAOD += metFilters


from MuonPhotonMatch_cff import muonPhotonMatch, muonPhotonMatchMiniAOD

leptons = cms.EDProducer('CustoTnPLeptonProducer_miniAOD',
                              muon_src = cms.InputTag('slimmedMuons'),
                              muon_srcSecond = cms.InputTag('slimmedMuons'),
                              muon_cuts = cms.string(''),
                              muon_track_for_momentum = cms.string('TunePNew'),
                              muon_track_for_momentum_CSC = cms.string('Inner'),
                              muon_photon_match_src = cms.InputTag('muonPhotonMatchMiniAOD'),
                              electron_muon_veto_dR = cms.double(-1),
                              trigger_match_max_dR = cms.double(0.2),
                              trigger_summary = cms.InputTag('slimmedPatTrigger'), # to run on 2017 Data
                              bits = cms.InputTag("TriggerResults","","HLT"),#data
                              prescales = cms.InputTag("patTrigger"),
                              )

leptonsAOD = cms.EDProducer('CustoTnPLeptonProducer_AOD',
                                muon_src = cms.InputTag('slimmedMuonsCustom'),
                                muon_srcSecond = cms.InputTag('slimmedMuonsCustom'),
                                trigger_summary = cms.InputTag('selectedPatTriggerCustom'),
                                muon_photon_match_src = cms.InputTag('muonPhotonMatchAOD'),
                                muon_cuts = cms.string(''),
                                muon_track_for_momentum = cms.string('TunePNew'),
                                muon_track_for_momentum_CSC = cms.string('Inner'),
                                electron_muon_veto_dR = cms.double(-1),
                                trigger_match_max_dR = cms.double(0.2),
                                bits = cms.InputTag("TriggerResults","","HLT"),
                                prescales = cms.InputTag("patTrigger"),
                                )

