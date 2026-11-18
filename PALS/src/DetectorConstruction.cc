//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the B3::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4NistManager.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"

namespace B3
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();

  G4bool isotopes = false;

  G4Element* O = man->FindOrBuildElement("O", isotopes);
  G4Element* Si = man->FindOrBuildElement("Si", isotopes);
  G4Element* Lu = man->FindOrBuildElement("Lu", isotopes);

  auto LSO = new G4Material("Lu2SiO5", 7.4 * g / cm3, 3);
  LSO->AddElement(Lu, 2);
  LSO->AddElement(Si, 1);
  LSO->AddElement(O, 5);

  //New materials start here

  //CeBr3
  G4Element* Ce = man->FindOrBuildElement("Ce", isotopes);
  G4Element* Br = man->FindOrBuildElement("Br", isotopes);

  auto CeBr3 = new G4Material("CeBr3", 5.10 * g / cm3, 2);
  CeBr3->AddElement(Ce, 1);
  CeBr3->AddElement(Br, 3);

  //USE UNDOPED NaI, CAN REVISIT LATER
  //NaI(Tl) (normal NaI does not have dopant)
  //G4Element* Na = man->FindOrBuildElement("Na", isotopes);
  //G4Element* I = man->FindOrBuildElement("I", isotopes);
  //G4Element* Tl = man->FindOrBuildElement("Tl", isotopes);

  //auto NaITl = new G4Material("NaI(Tl)", 3.67 * g / cm3, 2);
  //NaITl->AddElement(Na, 1);
  //NaITl->AddElement(I, 3);
  //NaITl->AddElement(Br, 3);

  //BC408
  G4Element* H = man->FindOrBuildElement("H", isotopes);
  G4Element* C = man->FindOrBuildElement("C", isotopes);

  auto BC408 = new G4Material("BC408", 1.023 * g / cm3, 2);

  BC408->AddElement(C, 0.915206); //use mass fractions
  BC408->AddElement(H, 0.084794);

  auto pterph = new G4Material("pterph", 1.23 * g / cm3, 2);

  pterph->AddElement(C, 0.9387); //use mass fractions
  pterph->AddElement(H, 0.0613);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Gamma detector Parameters
  //
  G4double cryst_dX = 6 * cm, cryst_dY = 6 * cm, cryst_dZ = 3 * cm;
  G4int nb_cryst = 32;
  G4int nb_rings = 9;
  //
  G4double dPhi = twopi / nb_cryst, half_dPhi = 0.5 * dPhi;
  G4double cosdPhi = std::cos(half_dPhi);
  G4double tandPhi = std::tan(half_dPhi);
  //
  G4double ring_R1 = 0.5 * cryst_dY / tandPhi;
  G4double ring_R2 = (ring_R1 + cryst_dZ) / cosdPhi;
  //
  G4double detector_dZ = nb_rings * cryst_dX;
  //
  G4NistManager* nist = G4NistManager::Instance(); //Instantiates the NIST finder
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* cryst_mat = nist->FindOrBuildMaterial("Lu2SiO5");

  //Get materials for our project
  G4Material* lead_mat = nist->FindOrBuildMaterial("G4_Pb"); //lead
  //G4Material* mount = nist->FindOrBuildMaterial(); //plastic mounting piece (TODO get material) 
  G4Material* cryst_analysis = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE"); //NaI (don't need to include Tl)
  G4Material* cryst_scatt = nist->FindOrBuildMaterial("CeBr3"); //NaI (don't need to include Tl)


  //
  // World
  //
  G4double world_sizeXY = 100 * cm; //100 cm world size (length of side)
  G4double world_sizeZ = 100 * cm; //100 cm world size (length of side)

  auto solidWorld =
    new G4Box("World",  // its name
              0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                        default_mat,  // its material
                                        "World");  // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
                                     G4ThreeVector(),  // at (0,0,0)
                                     logicWorld,  // its logical volume
                                     "World",  // its name
                                     nullptr,  // its mother  volume
                                     false,  // no boolean operation
                                     0,  // copy number
                                     fCheckOverlaps);  // checking overlaps

auto leadVis = new G4VisAttributes(G4Colour(0.8, 0.2, 0.2)); // Red color

// ==========================================
  // TYPE 1: Brick Type 1 (e.g., Side Walls along X)
  // ==========================================
  G4double b1_dX = 1.0 * 2.54 * cm;   // Thickness
  G4double b1_dY = 8.0 * 2.54 * cm;  // Height along Y (vertical)
  G4double b1_dZ = 6.75 * 2.54 * cm;  // Length along Z

  auto solidBrickType1 = new G4Box("BrickType1_Solid", 0.5 * b1_dX, 0.5 * b1_dY, 0.5 * b1_dZ);
  auto logicBrickType1 = new G4LogicalVolume(solidBrickType1, lead_mat, "BrickType1_LV");
  logicBrickType1->SetVisAttributes(leadVis);

  double lead1_xoff = 4.0 * 0.5 * 2.54 * cm + b1_dX*0.5; //Distance of thick lead brick from center

  // Place Type 1 Bricks on the -X and +X sides
  new G4PVPlacement(nullptr, G4ThreeVector(-lead1_xoff, 0.0, 0.0), logicBrickType1, "Brick_MinusX", logicWorld, false, 0, fCheckOverlaps);
  new G4PVPlacement(nullptr, G4ThreeVector( lead1_xoff, 0.0, 0.0), logicBrickType1, "Brick_PlusX",  logicWorld, false, 1, fCheckOverlaps);

  // ==========================================
  // TYPE 2: Brick Type 2 (e.g., Front/Back Walls along Z)
  // ==========================================
  G4double b2_dX = 4.0 * 2.54 * cm;  // Width
  G4double b2_dY = 8.0 * 2.54 * cm;  // Height along Y (vertical)
  G4double b2_dZ = 2.0 * 2.54 * cm;   // Thickness

  auto solidBrickType2 = new G4Box("BrickType2_Solid", 0.5 * b2_dX, 0.5 * b2_dY, 0.5 * b2_dZ);
  auto logicBrickType2 = new G4LogicalVolume(solidBrickType2, lead_mat, "BrickType2_LV");
  logicBrickType2->SetVisAttributes(leadVis);


  // Define the Cylindrical "Cutter" for the hole
  // Arguments for G4Tubs: name, innerRadius, outerRadius, halfZ, startAngle, spanAngle
  G4double holeRadius = 1.0 * cm; 
  G4double holeHalfZ  = (b2_dZ / 2.) * 1.5; // Make this slightly thicker than the brick's half-thickness (3.75 cm) to avoid surface artifacts
  
  auto solidCylHole = new G4Tubs("CylHole", 
                                 0.0,            // Inner radius (solid cylinder)
                                 holeRadius,     // Outer radius
                                 holeHalfZ,      // Half-length along Z
                                 0.0,            // Starting phi angle
                                 360.0 * deg);   // Full circle span

  // Rotate the Cylinder Cutter so its length aligns with the Z-axis
  // By default, G4Tubs extrudes along the Z-axis, which matches our brick thickness direction.
  G4RotationMatrix* noRot = nullptr;
  G4ThreeVector holeOffset(0.0, 0.0, 0.0);

  // Perform the Subtraction: Brick minus Cylinder
  auto solidBrickWithHole = new G4SubtractionSolid("BrickZ_WithCylHole", 
                                                    solidBrickType2, 
                                                    solidCylHole, 
                                                    noRot, 
                                                    holeOffset);

  // Wrap into a Logical Volume and Place It
  auto logicBrickWithHole = new G4LogicalVolume(solidBrickWithHole, lead_mat, "BrickZ_CylHole_LV");
  logicBrickWithHole->SetVisAttributes(leadVis);

  double lead2_zoff = 3.5*cm + b2_dZ*0.5; //Distance of thick lead brick from center

  // Place the perforated bricks along the Z-axis
  new G4PVPlacement(nullptr,G4ThreeVector(0.0, 0.0, -lead2_zoff),logicBrickWithHole,"LeadBrick_Z_CylHole", 
                    logicWorld, 
                    false, 
                    2, 
                    fCheckOverlaps);

  new G4PVPlacement(nullptr, 
                    G4ThreeVector(0.0, 0.0, lead2_zoff), 
                    logicBrickWithHole, 
                    "LeadBrick_Z_CylHole", 
                    logicWorld, 
                    false, 
                    2, 
                    fCheckOverlaps);


  //
  // ring
  //
  auto solidRing = new G4Tubs("Ring", ring_R1, ring_R2, 0.5 * cryst_dX, 0., twopi);

  auto logicRing = new G4LogicalVolume(solidRing,  // its solid
                                       default_mat,  // its material
                                       "Ring");  // its name

  //
  // define crystal
  //
  G4double gap = 0.5 * mm;  // a gap for wrapping
  G4double dX = cryst_dX - gap, dY = cryst_dY - gap;
  auto solidCryst = new G4Box("crystal", dX / 2, dY / 2, cryst_dZ / 2);

  auto logicCryst = new G4LogicalVolume(solidCryst,  // its solid
                                        cryst_mat,  // its material
                                        "CrystalLV");  // its name

/*  // place crystals within a ring
  //
  for (G4int icrys = 0; icrys < nb_cryst; icrys++) {
    G4double phi = icrys * dPhi;
    G4RotationMatrix rotm = G4RotationMatrix();
    rotm.rotateY(90 * deg);
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi), std::sin(phi), 0.);
    G4ThreeVector position = (ring_R1 + 0.5 * cryst_dZ) * uz;
    G4Transform3D transform = G4Transform3D(rotm, position);

    new G4PVPlacement(transform,  // rotation,position
                      logicCryst,  // its logical volume
                      "crystal",  // its name
                      logicRing,  // its mother  volume
                      false,  // no boolean operation
                      icrys,  // copy number
                      fCheckOverlaps);  // checking overlaps
  }

  //
  // full detector
  //
  auto solidDetector = new G4Tubs("Detector", ring_R1, ring_R2, 0.5 * detector_dZ, 0., twopi);

  auto logicDetector = new G4LogicalVolume(solidDetector,  // its solid
                                           default_mat,  // its material
                                           "Detector");  // its name

  //
  // place rings within detector
  //
  G4double OG = -0.5 * (detector_dZ + cryst_dX);
  for (G4int iring = 0; iring < nb_rings; iring++) {
    OG += cryst_dX;
    new G4PVPlacement(nullptr,  // no rotation
                      G4ThreeVector(0, 0, OG),  // position
                      logicRing,  // its logical volume
                      "ring",  // its name
                      logicDetector,  // its mother  volume
                      false,  // no boolean operation
                      iring,  // copy number
                      fCheckOverlaps);  // checking overlaps
  }

  //
  // place detector in world
  //
  new G4PVPlacement(nullptr,  // no rotation
                    G4ThreeVector(),  // at (0,0,0)
                    logicDetector,  // its logical volume
                    "Detector",  // its name
                    logicWorld,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    fCheckOverlaps);  // checking overlaps*/

  //
  // patient
  //
  G4double patient_radius = 8 * cm;
  G4double patient_dZ = 10 * cm;
  G4Material* patient_mat = nist->FindOrBuildMaterial("G4_BRAIN_ICRP");

  auto solidPatient = new G4Tubs("Patient", 0., patient_radius, 0.5 * patient_dZ, 0., twopi);

  auto logicPatient = new G4LogicalVolume(solidPatient,  // its solid
                                          patient_mat,  // its material
                                          "PatientLV");  // its name

  //
  // place patient in world
  //
  new G4PVPlacement(nullptr,  // no rotation
                    G4ThreeVector(),  // at (0,0,0)
                    logicPatient,  // its logical volume
                    "Patient",  // its name
                    logicWorld,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    fCheckOverlaps);  // checking overlaps

  // Visualization attributes
  //
  //logicRing->SetVisAttributes(G4VisAttributes::GetInvisible());
  //logicDetector->SetVisAttributes(G4VisAttributes::GetInvisible());

  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible()); //make world invisible

  // Create visibility attributes for the world
  auto worldVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); // Gray wireframe color
  worldVis->SetVisibility(true);
  worldVis->SetForceWireframe(true); // Forces it to render as lines only, never solid
  
  // Apply it to your world logical volume
  logicWorld->SetVisAttributes(worldVis);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  // always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // declare crystal as a MultiFunctionalDetector scorer
  //
  auto cryst = new G4MultiFunctionalDetector("crystal");
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
  SetSensitiveDetector("CrystalLV", cryst);

  // declare patient as a MultiFunctionalDetector scorer
  //
  auto patient = new G4MultiFunctionalDetector("patient");
  G4SDManager::GetSDMpointer()->AddNewDetector(patient);
  G4VPrimitiveScorer* primitiv2 = new G4PSDoseDeposit("dose");
  patient->RegisterPrimitive(primitiv2);
  SetSensitiveDetector("PatientLV", patient);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B3
