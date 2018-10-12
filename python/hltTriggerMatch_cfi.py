import FWCore.ParameterSet.Config as cms

def make_string_cut_for_trigger_matching( list_path_names, list_filters_pt ):
  cut = ''
  if len(list_path_names) != len(list_filters_pt):
    print 'len(list_path_names) != len(list_filters_pt) -> return ', cut
    return cut
  for i, f in enumerate(list_path_names):
    if f != list_path_names[-1]:
      cut += 'userFloat("%s_TriggerMatchPt")>=%i || ' % (list_path_names[i], list_filters_pt[i])
    else:
      cut += 'userFloat("%s_TriggerMatchPt")>=%i ' % (list_path_names[i], list_filters_pt[i])
  return cut


trigger_path_names_2016 = [
  'Mu50',
  'TkMu50'
]
trigger_filters_2016 = [
  'hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q',
  'hltL3fL1sMu25f0TkFiltered50Q'
]
trigger_filters_pt_2016 = [
  50,
  50
]

trigger_match_2018 = make_string_cut_for_trigger_matching( trigger_path_names_2016, trigger_filters_pt_2016 )
trigger_match_2016 = trigger_match_2018



# trigger_pt_threshold = 50
# offline_pt_threshold = 53 #?
# trigger_paths = ['HLT_Mu50_v%i' % i for i in (6, 7, 8, 9, 10, 11)]
# trigger_match = 'userFloat("TriggerMatchPt") > %(trigger_pt_threshold)i ' % locals()
#trigger_match = '1>0'

# overall_prescale = 350
# prescaled_trigger_pt_threshold = 24
# prescaled_offline_pt_threshold = 27
# prescaled_trigger_paths = ['HLT_Mu27_v%i' % i for i in (6, 7, 8, 9, 10, 11)]
# prescaled_trigger_match = trigger_match.replace('Trigger', 'prescaledTrigger').replace('%i' % trigger_pt_threshold, '%i' % prescaled_trigger_pt_threshold)








#--- Not using
# muonTriggerMatchHLTMuons = cms.EDProducer('PATTriggerMatcherDRLessByR',
#                                           src = cms.InputTag( 'cleanPatMuons' ),
#                                           matched = cms.InputTag( 'patTrigger' ),
#                                           matchedCuts = cms.string('type("TriggerMuon") && (path("HLT_Mu9*",1,0) || path("HLT_Mu15*",1,0) || path("HLT_Mu24_v*",1,0)|| path("HLT_Mu24*",1,0) || path("HLT_Mu30*",1,0) || path("HLT_Mu40*",1,0) || path("HLT_Mu45*",1,0) || path("HLT_Mu50*",1,0))'),
#                                           maxDPtRel   = cms.double( 1. ), 
#                                           maxDeltaR   = cms.double( 0.2 ),
#                                           resolveAmbiguities    = cms.bool( True ),
#                                           resolveByMatchQuality = cms.bool( True )
#                                           )

# muonTriggerMatchHLTMuonsMiniAOD = cms.EDProducer('PATTriggerMatcherDRLessByR',
#                                           src = cms.InputTag( 'slimmedMuons' ),
#                                           matched = cms.InputTag( 'patTrigger' ),
#                                           matchedCuts = cms.string('type("TriggerMuon") && (path("HLT_Mu9*",1,0) || path("HLT_Mu15*",1,0) || path("HLT_Mu24_v*",1,0)|| path("HLT_Mu24*",1,0) || path("HLT_Mu30*",1,0) || path("HLT_Mu40*",1,0) || path("HLT_Mu45*",1,0) || path("HLT_Mu50*",1,0))'),
#                                           maxDPtRel   = cms.double( 1. ), 
#                                           maxDeltaR   = cms.double( 0.2 ),
#                                           resolveAmbiguities    = cms.bool( True ),
#                                           resolveByMatchQuality = cms.bool( True )
# )
