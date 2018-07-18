import FWCore.ParameterSet.Config as cms

# The below is mostly copied from
# https://twiki.cern.ch/twiki/bin/view/CMS/Collisions2010Recipes ; the
# user is responsible for checking that what's used is up-to-date and
# appropriate to their analysis.

primaryVertex = cms.EDFilter('GoodVertexFilter',
                                   vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                   minimumNDOF = cms.uint32(4),
                                   maxAbsZ = cms.double(24),
                                   maxd0 = cms.double(2)
                                   )


from METFilter_miniAOD_cfi import defaultSelector

primaryVertexMiniAOD = defaultSelector.clone()
primaryVertexMiniAOD.flag = "Flag_goodVertices"

#beamHaloMiniAOD = defaultSelector.clone()
#beamHaloMiniAOD.flag = "Flag_globalTightHalo2016Filter"

#hbheNoiseMiniAOD = defaultSelector.clone()
#hbheNoiseMiniAOD.flag = "Flag_HBHENoiseFilter"

#hbheIsoNoiseMiniAOD = defaultSelector.clone()
#hbheIsoNoiseMiniAOD.flag = "Flag_HBHENoiseIsoFilter"

#ecalDeadCellPrimitiveMiniAOD = defaultSelector.clone()
#ecalDeadCellPrimitiveMiniAOD.flag = "Flag_EcalDeadCellTriggerPrimitiveFilter"

#eeBadScMiniAOD = defaultSelector.clone()
#eeBadScMiniAOD.flag = "Flag_eeBadScFilter"

#metFilters = [beamHaloMiniAOD, hbheNoiseMiniAOD, hbheIsoNoiseMiniAOD, ecalDeadCellPrimitiveMiniAOD, eeBadScMiniAOD]

'''
def addNewFilters(process,filterList):

	process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
	process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
	process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

	process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
	process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
	process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

	filterList.append(process.BadPFMuonFilter)
	filterList.append(process.BadChargedCandidateFilter)

	return filterList
'''

