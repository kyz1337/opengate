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
                    const std::vector<std::vector<int>> &,
                    const std::vector<std::vector<G4double>> &,
                    bool,
                    G4double>(),
           py::arg("xs_scaling"),
           py::arg("desired_channel"),
           py::arg("energy_ranges")            = std::vector<std::vector<G4double>>{},
           py::arg("exclusive")                = false,
           py::arg("projectile_cut_fraction")  = 0.0,
           "Construct the physics constructor.\n\n"
           "xs_scaling               -- factor applied to the total alphaInelastic cross section.\n"
           "desired_channel          -- list of [Z, A] pairs that must appear in the\n"
           "                           reaction products.\n"
           "                           Example He3+n+X: [[2,3],[0,1]]\n"
           "energy_ranges            -- list of [Z, A, E_min, E_max] constraints (MeV);\n"
           "                           use -1 for no bounds.\n"
           "exclusive                -- if set to True, use exclusive channel definition (exact match)"
           "projectile_cut_fraction  -- if > 0, automatically require composite\n"
           "                           fragments (A>=2) in desired_channel to have\n"
           "                           KE >= fraction * A * KE_beam/4, separating\n"
           "                           projectile from target fragments dynamically.\n"
           "                           0 = disabled (default).");
}
