//
//  Implementation of modular detector subroutines for GEANT4
//

#include "DetectorModules.hh"

#include <algorithm>

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4String.hh"


namespace {

//Specifically for cylinder-type detectors with an Al can
G4LogicalVolume* BuildCylAssemblyLV(const G4String& prefix,
                                    G4Material* crystalMat,
                                    G4Material* alMat,
                                    G4Material* glassMat,
                                    G4Material* envMat,
                                    G4double cylRadius,
                                    G4double halfZ,
                                    G4double alThick,
                                    G4double winThick,
                                    G4bool checkOverlaps,
                                    const G4Colour& crystalColor) {

  //Get half-zs of can and glass window
  const G4double capHalfZ = 0.5 * alThick;
  const G4double winHalfZ = 0.5 * winThick;

  // Mother assembly volume enclosing crystal, front-capped Al can, and rear glass window
  // 0.01 refers to the snugness of the parts
  const G4double asmR       = cylRadius + alThick + 0.01 * mm;
  const G4double maxExtentZ = std::max(halfZ + alThick, halfZ + winThick);
  const G4double asmHalfZ   = maxExtentZ + 0.01 * mm;

  //Generates solid and logical volumes for the assembly, contains all parts
  auto solidAsm = new G4Tubs(prefix + "_Assembly_Solid",
                             0.0,
                             asmR,
                             asmHalfZ,
                             0.0,
                             360.0 * deg);
  auto logicAsm = new G4LogicalVolume(solidAsm, envMat, prefix + "_Assembly_LV");

  //Colors the assemblies
  auto asmVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.15));
  asmVis->SetForceWireframe(true);
  logicAsm->SetVisAttributes(asmVis);

  //Scintillation crystal centered at local (0, 0, 0)
  auto solidCrystal = new G4Tubs(prefix + "_Det_Solid",
                                 0.0,
                                 cylRadius,
                                 halfZ,
                                 0.0,
                                 360.0 * deg);
  auto logicCrystal = new G4LogicalVolume(solidCrystal, crystalMat, prefix + "_Det_LV");

  //Color the crystal
  auto crystalVis = new G4VisAttributes(crystalColor);
  crystalVis->SetForceSolid(true);
  logicCrystal->SetVisAttributes(crystalVis);

  //Places the crystal in the volume
  new G4PVPlacement(nullptr,
                    G4ThreeVector(0.0, 0.0, 0.0),
                    logicCrystal,
                    prefix + "_Crystal",
                    logicAsm,
                    false,
                    0,
                    checkOverlaps);

  //Al coloring
  auto alVis = new G4VisAttributes(G4Colour(0.75, 0.75, 0.8, 0.5));
  alVis->SetForceSolid(true);

  //Generates Aluminum can that fits around the crystal (not including the face)
  auto solidAlWall = new G4Tubs(prefix + "_AlWall_Solid",
                                cylRadius,
                                cylRadius + alThick,
                                halfZ,
                                0.0,
                                360.0 * deg);
  auto logicAlWall = new G4LogicalVolume(solidAlWall, alMat, prefix + "_AlWall_LV");
  logicAlWall->SetVisAttributes(alVis);

  //Places the Al can
  new G4PVPlacement(nullptr,
                    G4ThreeVector(0.0, 0.0, 0.0),
                    logicAlWall,
                    prefix + "_AlWall",
                    logicAsm,
                    false,
                    0,
                    checkOverlaps);

  //Caps the front of the crystal with an Al disk
  const G4double frontCapZ = halfZ + capHalfZ;
  auto solidAlCap = new G4Tubs(prefix + "_AlFrontCap_Solid",
                               0.0,
                               cylRadius + alThick,
                               capHalfZ,
                               0.0,
                               360.0 * deg);
  auto logicAlCap = new G4LogicalVolume(solidAlCap, alMat, prefix + "_AlFrontCap_LV");
  logicAlCap->SetVisAttributes(alVis);

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0.0, 0.0, frontCapZ),
                    logicAlCap,
                    prefix + "_AlFrontCap",
                    logicAsm,
                    false,
                    0,
                    checkOverlaps);

  //PMT glass window on the rear (-Z) face
  const G4double winLocalZ = -halfZ - winHalfZ;
  auto solidWin = new G4Tubs(prefix + "_GlassWin_Solid",
                             0.0,
                             cylRadius,
                             winHalfZ,
                             0.0,
                             360.0 * deg);
  auto logicWin = new G4LogicalVolume(solidWin, glassMat, prefix + "_GlassWin_LV");

  auto winVis = new G4VisAttributes(G4Colour(0.4, 0.8, 1.0, 0.45));
  winVis->SetForceSolid(true);
  logicWin->SetVisAttributes(winVis);

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0.0, 0.0, winLocalZ),
                    logicWin,
                    prefix + "_GlassWin",
                    logicAsm,
                    false,
                    0,
                    checkOverlaps);
  //Return logical volume of the whole assembly
  return logicAsm;
}

} //end namespace

//Now call the individual assemblies
namespace DetectorModules {

  //Standard detector assembly
  G4LogicalVolume* BuildAnalyzerAssemblyLV(const G4String& prefix,
                                            G4Material* crystalMat,
                                            G4Material* alMat,
                                            G4Material* glassMat,
                                            G4Material* envMat,
                                            const DetCylinder& cfg) {
  return BuildCylAssemblyLV(prefix,
                            crystalMat,
                            alMat,
                            glassMat,
                            envMat,
                            cfg.cylRadius,
                            cfg.halfZ,
                            cfg.alThick,
                            cfg.winThick,
                            cfg.checkOverlaps,
                            cfg.crystalColor);
  }


} //end DetectorModules
