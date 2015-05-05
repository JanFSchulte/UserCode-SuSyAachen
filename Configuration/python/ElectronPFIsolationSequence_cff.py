import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Isolation.tools_cfi import *

def load_electronPFiso_sequence(proc, seq_name, algo, coneR, src, src_charged_hadron='', src_neutral_hadron='', src_photon='', src_charged_pileup=''):

    doCH, doNH, doPh, doPU = False, False, False, False
    if src_charged_hadron != '': doCH = True
    if src_neutral_hadron != '': doNH = True
    if src_photon         != '': doPh = True
    if src_charged_pileup != '': doPU = True

    iso_seq = cms.Sequence()

    if doCH:
        setattr(proc, 'elePFIsoDepositCH'+algo, isoDepositReplace(src, src_charged_hadron))
        iso_seq += getattr(proc, 'elePFIsoDepositCH'+algo)

    if doNH:
        setattr(proc, 'elePFIsoDepositNH'+algo, isoDepositReplace(src, src_neutral_hadron))
        iso_seq += getattr(proc, 'elePFIsoDepositNH'+algo)

    if doPh:
        setattr(proc, 'elePFIsoDepositPh'+algo, isoDepositReplace(src, src_photon))
        iso_seq += getattr(proc, 'elePFIsoDepositPh'+algo)

    if doPU:
        setattr(proc, 'elePFIsoDepositPU'+algo, isoDepositReplace(src, src_charged_pileup))
        iso_seq += getattr(proc, 'elePFIsoDepositPU'+algo)

    iso_vals_seq = cms.Sequence()

    if doCH:
        setattr(proc, 'elePFIsoValueCH'+algo,
          cms.EDProducer('CandIsolatorFromDeposits',
            deposits = cms.VPSet(
              cms.PSet(
                src = cms.InputTag('elePFIsoDepositCH'+algo),
                deltaR = cms.double(coneR),
                weight = cms.string('1'),
                vetos = cms.vstring('0.0001','Threshold(0.0)'),
                skipDefaultVeto = cms.bool(True),
                mode = cms.string('sum')
              )
            )
          )
        )
        iso_vals_seq += getattr(proc, 'elePFIsoValueCH'+algo)

    if doNH:
        setattr(proc, 'elePFIsoValueNH'+algo,
          cms.EDProducer('CandIsolatorFromDeposits',
            deposits = cms.VPSet(
              cms.PSet(
                src = cms.InputTag('elePFIsoDepositNH'+algo),
                deltaR = cms.double(coneR),
                weight = cms.string('1'),
                vetos = cms.vstring('0.01','Threshold(0.5)'),
                skipDefaultVeto = cms.bool(True),
                mode = cms.string('sum')
              )
            )
          )
        )
        iso_vals_seq += getattr(proc, 'elePFIsoValueNH'+algo)

    if doPh:
        setattr(proc, 'elePFIsoValuePh'+algo,
          cms.EDProducer('CandIsolatorFromDeposits',
            deposits = cms.VPSet(
              cms.PSet(
                src = cms.InputTag('elePFIsoDepositPh'+algo),
                deltaR = cms.double(coneR),
                weight = cms.string('1'),
                vetos = cms.vstring('0.01','Threshold(0.5)'),
                skipDefaultVeto = cms.bool(True),
                mode = cms.string('sum')
              )
            )
          )
        )
        iso_vals_seq += getattr(proc, 'elePFIsoValuePh'+algo)

    if doPU:
        setattr(proc, 'elePFIsoValuePU'+algo,
          cms.EDProducer('CandIsolatorFromDeposits',
            deposits = cms.VPSet(
              cms.PSet(
                src = cms.InputTag('elePFIsoDepositPU'+algo),
                deltaR = cms.double(coneR),
                weight = cms.string('1'),
                vetos = cms.vstring('0.01','Threshold(0.5)'),
                skipDefaultVeto = cms.bool(True),
                mode = cms.string('sum')
              )
            )
          )
        )
        iso_vals_seq += getattr(proc, 'elePFIsoValuePU'+algo)

    iso_seq *= iso_vals_seq

    setattr(proc, seq_name, iso_seq)
