//
//  Addon file for GEANT4 that handles the construction of modular detectors
//  Written by Johnathan Phillips
//  11/20/2026
//
//  Detector parts are contained within structs and build by subroutines
//
/// \file DetectorConstruction.hh
/// \brief Definition of the B3::DetectorConstruction class

#pragma once

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"

//Forward declarations
class G4LogicalVolume;
class G4Material;

namespace DetectorModules {

//generic struct for cylindrical detectors with Al can and PMT glass
struct DetCylinder {
  DetCylinder(G4double _cylRadius,G4double _halfZ,G4double _alThick,G4double _winThick, G4bool _checkOverlaps,G4Colour _crystalColor) {
    cylRadius = _cylRadius;
    halfZ = _halfZ;
    alThick = _alThick;
    winThick = _winThick;
    checkOverlaps = _checkOverlaps;
    crystalColor = _crystalColor;
  }
  G4double cylRadius; //Radius of cylinder
  G4double halfZ; //Distance from Z-center to one end
  G4double alThick; //Thickness of aluminum can
  G4double winThick; //thickness of phototube glass
  G4bool checkOverlaps = false;
  G4Colour crystalColor; //Color of the crystal
};
  
  //Standard teaching lab NaI detectors
  const DetCylinder NaIAnalyzerConfig{12.5 * mm, 12.5 * mm, 0.1  * mm, 8.0  * mm, false,G4Colour(0.45, 0.25, 0.0)};

  //Scionix 25 mm diameter CeBr3 detectors
  const DetCylinder CeBr3_25_AnalyzerConfig{12.5 * mm, 12.5 * mm, 0.4  * mm, 8.0  * mm, false,G4Colour(0.0, 0.85, 0.95)};

// Builds the 3-part cylinder assembly
G4LogicalVolume* BuildAnalyzerAssemblyLV(const G4String& prefix,
                                            G4Material* crystalMat,
                                            G4Material* alMat,
                                            G4Material* glassMat,
                                            G4Material* envMat,
                                            const DetCylinder& cfg);

}


/*//Standard NaI detectors in the teaching lab
struct NaIAnalyzerConfig {
  const G4double cylRadius     = 12.5 * mm; //Radius of cylinder
  const G4double halfZ         = 12.5 * mm; //Distance from Z-center to one end
  const G4double alThick       = 0.1  * mm; //Thickness of aluminum can
  const G4double winThick      = 8.0  * mm; //thickness of phototube glass
  G4bool   checkOverlaps = false;
};

//25 mm diameter CeBr3 detectors from Scionix
struct CeBr3_25_AnalyzerConfig {
  const G4double cylRadius     = 12.5 * mm; //Radius of cylinder
  const G4double halfZ         = 12.5 * mm; //Distance from Z-center to one end
  const G4double alThick       = 0.4  * mm; //Thickness of aluminum can
  const G4double winThick      = 8.0  * mm; //thickness of phototube glass
  G4bool   checkOverlaps = false;
};*/

