#pragma once

#include "globals.hh"

struct AnalyzerPlacement {
  G4int id;
  G4String module;
  G4String side;
  G4double distance;  // Scatterer center to outer front face (aluminum cap).
  G4double theta;
  G4double phi;
};
