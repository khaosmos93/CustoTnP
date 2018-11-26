import FWCore.ParameterSet.Config as cms

# A filter for post-tuple filtering on the goodData results as stored
# in a TriggerResults object instead of filtering at tuple-making
# time.
'''
import HLTrigger.HLTfilters.hltHighLevel_cfi
goodDataFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
goodDataFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', 'PAT')
goodDataFilter.HLTPaths = ["goodDataPrimaryVertexFilter"] # can set to just 'goodDataPrimaryVertexFilter', for example
#goodDataFilter.HLTPaths = ['goodDataMETFilter']
goodDataFilter.andOr = False # = AND
'''

from CustoTnP.Analyzer.goodData_cff import primaryVertexMiniAOD #, hltPhysicsDeclared, metFilters
goodDataFiltersMiniAOD = [primaryVertexMiniAOD]

## for full filtering, use:
#goodDataFiltersMiniAOD = [primaryVertexMiniAOD,hltPhysicsDeclared]
#goodDataFiltersMiniAOD += metFilters


from MuonPhotonMatch_cff import muonPhotonMatch, muonPhotonMatchMiniAOD
from CustoTnP.Analyzer.hltTriggerMatch_cfi import trigger_path_names, trigger_filters, trigger_filters_pt

leptons = cms.EDProducer('CustoTnPLeptonProducer',
                              muon_src = cms.InputTag('slimmedMuons'),
                              muon_srcSecond = cms.InputTag('slimmedMuons'),
                              muon_cuts = cms.string(''),
                              muon_track_for_momentum = cms.string('TunePNew'),
                              muon_track_for_momentum_CSC = cms.string('Inner'),
                              muon_photon_match_src = cms.InputTag('muonPhotonMatchMiniAOD'),
                              trigger_match_max_dR = cms.double(0.2),
                              trigger_summary = cms.InputTag('slimmedPatTrigger'), # to run on 2017 Data
                              bits = cms.InputTag("TriggerResults","","HLT"),
                              prescales = cms.InputTag("patTrigger"),
                              trigger_filters = cms.vstring(trigger_filters),
                              trigger_path_names = cms.vstring(trigger_path_names),
                              l1 = cms.InputTag("gmtStage2Digis", "Muon", "RECO"),
                              reco_muon_src = cms.InputTag(''),            # only for AOD
                              muonshower_src = cms.InputTag('', '', ''),   # only for AOD
)

leptonsAOD = cms.EDProducer('CustoTnPLeptonProducer',
                              muon_src = cms.InputTag('slimmedMuonsCustom'),
                              muon_srcSecond = cms.InputTag('slimmedMuonsCustom'),
                              muon_cuts = cms.string(''),
                              muon_track_for_momentum = cms.string('TunePNew'),
                              muon_track_for_momentum_CSC = cms.string('Inner'),
                              muon_photon_match_src = cms.InputTag('muonPhotonMatchAOD'),
                              trigger_match_max_dR = cms.double(0.2),
                              trigger_summary = cms.InputTag('selectedPatTriggerCustom'),
                              bits = cms.InputTag("TriggerResults","","HLT"),
                              prescales = cms.InputTag("patTrigger"),
                              trigger_filters = cms.vstring(trigger_filters),
                              trigger_path_names = cms.vstring(trigger_path_names),
                              l1 = cms.InputTag("gmtStage2Digis", "Muon", "RECO"),
                              reco_muon_src = cms.InputTag('muons'),
                              muonshower_src = cms.InputTag('muons', 'muonShowerInformation', 'RECO'),
)

