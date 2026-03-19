/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "G4VPhysicsConstructor.hh"
#include "GateChannelSelectiveWrapper.h"

void init_GateChannelSelectiveWrapperPhysics(py::module &m) {
  py::class_<GateChannelSelectiveWrapperPhysics, G4VPhysicsConstructor,
             std::unique_ptr<GateChannelSelectiveWrapperPhysics, py::nodelete>>(
      m, "GateChannelSelectiveWrapperPhysics")
      .def(py::init<G4double,
                    const std::string &,
                    const std::vector<std::vector<int>> &,
                    const std::vector<std::vector<G4double>> &,
                    bool>(),
           py::arg("xs_scaling"),
           py::arg("process_name"),
           py::arg("desired_channel"),
           py::arg("energy_ranges") = std::vector<std::vector<G4double>>{},
           py::arg("exclusive")     = false,
           "Construct the physics constructor.\n\n"
           "xs_scaling      -- factor applied to the hadronic inelastic cross section.\n"
           "process_name    -- Geant4 process name to wrap, e.g. 'alphaInelastic',\n"
           "                   'protonInelastic', 'He3Inelastic'. The owning particle\n"
           "                   is found automatically by searching all process managers.\n"
           "desired_channel -- list of [Z, A] pairs that must appear in the final state.\n"
           "                   Example He3+n+X: [[2,3],[0,1]]\n"
           "                   Note: for non-alpha projectiles, include the scattered\n"
           "                   primary in desired_channel if using exclusive=True.\n"
           "energy_ranges   -- list of [Z, A, E_min, E_max] constraints (MeV);\n"
           "                   use -1 for no bound on that side.\n"
           "exclusive       -- if True, final state must exactly match desired_channel.");
}
