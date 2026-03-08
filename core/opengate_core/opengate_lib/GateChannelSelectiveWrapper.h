/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateChannelSelectiveWrapper_h
#define GateChannelSelectiveWrapper_h

#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "G4ParticleChange.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4WrapperProcess.hh"

class GateChannelSelectiveWrapper : public G4WrapperProcess {
public:
  using ChannelSig   = std::vector<std::pair<int, int>>;
  using ChannelEntry = std::tuple<int, int, G4double>;
  using ChannelFull  = std::vector<ChannelEntry>;

  GateChannelSelectiveWrapper(
      G4VProcess *wrappedProcess, G4double xsScaling,
      const std::vector<std::vector<int>>    &desiredChannel,
      const std::vector<std::vector<G4double>> &energyRanges = {},
      bool exclusive = false,
      G4double projectileCutFraction = 0.0);
  ~GateChannelSelectiveWrapper() override;

  G4double PostStepGetPhysicalInteractionLength(const G4Track &track,
                                                G4double previousStepSize,
                                                G4ForceCondition *condition) override;

  G4VParticleChange *PostStepDoIt(const G4Track &track,
                                  const G4Step &step) override;

private:
  G4double   fXSScaling;
  ChannelSig fDesiredChannel;
  std::map<std::pair<int,int>, std::pair<G4double,G4double>> fEnergyRanges;
  bool     fExclusive;
  G4double fProjectileCutFraction;

  G4ParticleChange fNullChange;

  static ChannelFull BuildChannelFull(const G4VParticleChange *pc,
                                      const G4Track *track);

  bool IsDesiredChannel(const ChannelFull &full, G4double primaryKE) const;
};

class GateChannelSelectiveWrapperPhysics : public G4VPhysicsConstructor {
public:
  GateChannelSelectiveWrapperPhysics(
      G4double xsScaling,
      const std::vector<std::vector<int>>    &desiredChannel,
      const std::vector<std::vector<G4double>> &energyRanges = {},
      bool exclusive = false,
      G4double projectileCutFraction = 0.0);
  ~GateChannelSelectiveWrapperPhysics() override = default;

  void ConstructParticle() override {}
  void ConstructProcess() override;

private:
  G4double fXSScaling;
  std::vector<std::vector<int>>      fDesiredChannel;
  std::vector<std::vector<G4double>> fEnergyRanges;
  bool     fExclusive;
  G4double fProjectileCutFraction;
};

#endif
