/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateChannelSelectiveWrapper.h"

#include <algorithm>
#include <map>
#include <tuple>

#include "G4Alpha.hh"
#include "G4ForceCondition.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleChange.hh"
#include "G4VParticleChange.hh"
#include "Randomize.hh"


GateChannelSelectiveWrapper::GateChannelSelectiveWrapper(
    G4VProcess *wrappedProcess, G4double xsScaling,
    const std::vector<std::vector<int>>    &desiredChannel,
    const std::vector<std::vector<G4double>> &energyRanges,
    bool exclusive,
    G4double projectileCutFraction)
    : G4WrapperProcess("ChannelSelectiveWrapper_" +
                           wrappedProcess->GetProcessName(),
                       wrappedProcess->GetProcessType()),
      fXSScaling(xsScaling), fExclusive(exclusive),
      fProjectileCutFraction(projectileCutFraction) {
  RegisterProcess(wrappedProcess);

  for (const auto &pair : desiredChannel)
    fDesiredChannel.push_back({pair[0], pair[1]});
  std::sort(fDesiredChannel.begin(), fDesiredChannel.end());

  for (const auto &entry : energyRanges) {
    auto species = std::make_pair(static_cast<int>(entry[0]),
                                  static_cast<int>(entry[1]));
    fEnergyRanges[species] = {entry[2], entry[3]};
  }
}

GateChannelSelectiveWrapper::~GateChannelSelectiveWrapper() = default;

G4double GateChannelSelectiveWrapper::PostStepGetPhysicalInteractionLength(
    const G4Track &track, G4double previousStepSize,
    G4ForceCondition *condition) {
  if (fXSScaling < 1.0) {
    return G4WrapperProcess::PostStepGetPhysicalInteractionLength(
        track, previousStepSize, condition);
  }

  const G4double mfp =
      G4WrapperProcess::PostStepGetPhysicalInteractionLength(
          track, previousStepSize * fXSScaling, condition);

  return mfp / fXSScaling;
}

G4VParticleChange *
GateChannelSelectiveWrapper::PostStepDoIt(const G4Track &track,
                                           const G4Step &step) {
  G4VParticleChange *pc = G4WrapperProcess::PostStepDoIt(track, step);

  const ChannelFull full = BuildChannelFull(pc, &track);
  const G4double primaryKE = track.GetKineticEnergy();
  const bool isDesired = IsDesiredChannel(full, primaryKE);

  if (fXSScaling < 1.0) {
    if (!isDesired)
      return pc;
    if (G4UniformRand() <= fXSScaling)
      return pc;
    const G4int nSec = pc->GetNumberOfSecondaries();
    for (G4int i = 0; i < nSec; ++i) delete pc->GetSecondary(i);
    pc->Initialize(track);
    fNullChange.Initialize(track);
    return &fNullChange;
  }

  if (isDesired)
    return pc;

  if (G4UniformRand() * fXSScaling <= 1.0)
    return pc;

  const G4int nSec = pc->GetNumberOfSecondaries();
  for (G4int i = 0; i < nSec; ++i)
    delete pc->GetSecondary(i);
  pc->Initialize(track);
  fNullChange.Initialize(track);
  return &fNullChange;
}

GateChannelSelectiveWrapper::ChannelFull
GateChannelSelectiveWrapper::BuildChannelFull(const G4VParticleChange *pc,
                                               const G4Track *track) {
  ChannelFull full;

  auto isNuclear = [](const G4ParticleDefinition *def) -> bool {
    return (static_cast<int>(def->GetAtomicNumber()) > 0 ||
            static_cast<int>(def->GetAtomicMass()) > 0);
  };

  if (pc->GetTrackStatus() != fStopAndKill) {
    const auto *def = track->GetDefinition();
    if (isNuclear(def)) {
      const auto *pc4 = dynamic_cast<const G4ParticleChange *>(pc);
      const G4double primaryKE =
          pc4 ? pc4->GetEnergy() : track->GetKineticEnergy();
      full.push_back({static_cast<int>(def->GetAtomicNumber()),
                      static_cast<int>(def->GetAtomicMass()),
                      primaryKE});
    }
  }

  const G4int nSec = pc->GetNumberOfSecondaries();
  for (G4int i = 0; i < nSec; ++i) {
    const G4Track *sec = pc->GetSecondary(i);
    const auto *def = sec->GetDefinition();
    if (isNuclear(def))
      full.push_back({static_cast<int>(def->GetAtomicNumber()),
                      static_cast<int>(def->GetAtomicMass()),
                      sec->GetKineticEnergy()});
  }

  return full;
}

bool GateChannelSelectiveWrapper::IsDesiredChannel(
    const ChannelFull &full, G4double primaryKE) const {

  std::map<std::pair<int,int>, int> available;
  ChannelSig actualSig;
  for (const auto &[z, a, ke] : full) {
    auto sp = std::make_pair(z, a);
    available[sp]++;
    actualSig.push_back(sp);
  }
  std::sort(actualSig.begin(), actualSig.end());

  if (fExclusive) {
    if (actualSig != fDesiredChannel) return false;
  } else {
    std::map<std::pair<int,int>, int> required;
    for (const auto &r : fDesiredChannel) required[r]++;
    for (const auto &[species, count] : required) {
      if (available[species] < count) return false;
    }
  }

  for (const auto &[species, range] : fEnergyRanges) {
    int req_count = 0;
    for (const auto &r : fDesiredChannel) {
      if (r == species) req_count++;
    }
    if (req_count == 0) continue;

    const G4double lo = range.first;
    const G4double hi = range.second;

    int in_range = 0;
    for (const auto &[z, a, ke] : full) {
      if (std::make_pair(z, a) != species) continue;
      const bool ok = ((lo < 0.0) || (ke >= lo)) && ((hi < 0.0) || (ke <= hi));
      if (ok) in_range++;
    }
    if (in_range < req_count) return false;
  }

  if (fProjectileCutFraction > 0.0 && primaryKE > 0.0) {
    const G4double ke_per_u_beam = primaryKE / 4.0;

    std::map<std::pair<int,int>, int> required_proj;
    for (const auto &r : fDesiredChannel) {
      if (r.second >= 2) required_proj[r]++;
    }

    for (const auto &[species, req_count] : required_proj) {
      const G4double ke_thresh =
          fProjectileCutFraction * species.second * ke_per_u_beam;
      int n_above = 0;
      for (const auto &[z, a, ke] : full) {
        if (std::make_pair(z, a) != species) continue;
        if (ke >= ke_thresh) n_above++;
      }
      if (n_above < req_count) return false;
    }
  }

  return true;
}

GateChannelSelectiveWrapperPhysics::GateChannelSelectiveWrapperPhysics(
    G4double xsScaling,
    const std::vector<std::vector<int>>    &desiredChannel,
    const std::vector<std::vector<G4double>> &energyRanges,
    bool exclusive,
    G4double projectileCutFraction)
    : G4VPhysicsConstructor("GateChannelSelectiveWrapperPhysics"),
      fXSScaling(xsScaling), fDesiredChannel(desiredChannel),
      fEnergyRanges(energyRanges), fExclusive(exclusive),
      fProjectileCutFraction(projectileCutFraction) {}

void GateChannelSelectiveWrapperPhysics::ConstructProcess() {
  G4ParticleDefinition *alpha =
      G4ParticleTable::GetParticleTable()->FindParticle("alpha");
  G4ProcessManager *pm = alpha->GetProcessManager();

  G4VProcess *target = nullptr;
  G4ProcessVector *pvec = pm->GetProcessList();
  for (std::size_t i = 0; i < static_cast<std::size_t>(pvec->size()); ++i) {
    G4VProcess *p = (*pvec)[i];
    if (p && p->GetProcessName() == "alphaInelastic") {
      target = p;
      break;
    }
  }

  auto *wrapper = new GateChannelSelectiveWrapper(target, fXSScaling,
                                                   fDesiredChannel,
                                                   fEnergyRanges, fExclusive,
                                                   fProjectileCutFraction);
  pm->RemoveProcess(target);
  pm->AddDiscreteProcess(wrapper);
}
