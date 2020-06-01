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
//
/// \file B3DetectorConstruction.cc
/// \brief Implementation of the B3DetectorConstruction class

#include "B3DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::B3DetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true)
{
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::~B3DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();

  G4bool isotopes = false;

  G4Element*  O = man->FindOrBuildElement("O" , isotopes);
  G4Element* Si = man->FindOrBuildElement("Si", isotopes);
  G4Element* Lu = man->FindOrBuildElement("Lu", isotopes);

  G4Material* LSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
  LSO->AddElement(Lu, 2);
  LSO->AddElement(Si, 1);
  LSO->AddElement(O , 5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B3DetectorConstruction::Construct()
{
  // Gamma detector Parameters
  //
  G4double cryst_dX = 6*cm, cryst_dY = 6*cm, cryst_dZ = 3*cm;
  G4int nb_cryst = 32;
  G4int nb_rings = 9;
  //
  G4double dPhi = twopi/nb_cryst, half_dPhi = 0.5*dPhi;
  G4double cosdPhi = std::cos(half_dPhi);
  G4double tandPhi = std::tan(half_dPhi);
  //
  G4double ring_R1 = 0.5*cryst_dY/tandPhi;
  G4double ring_R2 = (ring_R1+cryst_dZ)/cosdPhi;
  //
  G4double detector_dZ = nb_rings*cryst_dX;
  //
  G4NistManager* nist = G4NistManager::Instance();

  G4Material* Vacuum = new G4Material("Galactic", 1., 1.01*g/mole,
                       1.e-25*g/cm3,kStateGas,2.73*kelvin,3.e-18*pascal);

  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* cryst_mat   = nist->FindOrBuildMaterial("G4_Ge");

  //
  // World
  //
  G4double world_sizeXY = 2.4*ring_R2;
  G4double world_sizeZ  = 1.2*detector_dZ;

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ); //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        default_mat,         //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps



  //
  // define crystal
  //
  G4double gap = 0.5*mm;        //a gap for wrapping
  G4double dX = cryst_dX - gap, dY = cryst_dY - gap;

  //user define params
  G4double distance_to_bewindow=50.*mm;
  G4double gap_bewindow_deadlayer=5.*mm;

  G4double bewindow_radius = 42*mm;
  G4double bewindow_dZ = 0.6*mm;
  G4double bewindow_Z = distance_to_bewindow;
  G4double deadlayer_radius = 42*mm;
  G4double deadlayer_dZ = 0.003*mm;
  G4double deadlayer_Z=bewindow_Z+bewindow_dZ/2+gap_bewindow_deadlayer+deadlayer_dZ/2;
  G4double hpge_radius=40.*mm;
  G4double hpge_dZ=30.*mm;
  G4double hpge_Z=deadlayer_Z+deadlayer_dZ/2+hpge_dZ/2;

  G4Tubs* solidCryst =
    new G4Tubs("crystal", 0, hpge_radius, hpge_dZ/2, 0., twopi);

  //G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, cryst_dZ/2);

  G4LogicalVolume* logicCryst =
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name

  //
  // place crystal in world
  //
  new G4PVPlacement(0,                     //no rotation
                    G4ThreeVector(0,0,hpge_Z), //position
                    logicCryst,             //its logical volume
                    "crystal",                //its name
                    logicWorld,         //its mother  volume
                    false,                 //no boolean operation
                    0,                 //copy number
                    fCheckOverlaps);       // checking overlaps

  //
  // be windows
  //
  G4Material* patient_mat = nist->FindOrBuildMaterial("G4_Be");

  G4Tubs* solidPatient =
    new G4Tubs("Patient", 0., bewindow_radius, 0.5*bewindow_dZ, 0., twopi);

  G4LogicalVolume* logicPatient =
    new G4LogicalVolume(solidPatient,        //its solid
                        patient_mat,         //its material
                        "PatientLV");        //its name

  //
  // place patient in world
  //
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,bewindow_Z),         //at (0,0,0)
                    logicPatient,            //its logical volume
                    "Patient",               //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps

  //
  // deadlayer windows
  //


  G4Tubs* solidDeadLayer =
    new G4Tubs("DeadLayer", 0., deadlayer_radius, 0.5*deadlayer_dZ, 0., twopi);

  G4LogicalVolume* logicDeadLayer =
    new G4LogicalVolume(solidDeadLayer,        //its solid
                        cryst_mat,         //its material
                        "DeadLayerLV");        //its name

  //
  // place patient in world
  //
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,deadlayer_Z),         //at (0,0,0)
                    logicDeadLayer,            //its logical volume
                    "DeadLayer",               //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps

  // Visualization attributes
  //
  G4VisAttributes* vishpge=new G4VisAttributes(G4Colour(1,0,0));
  vishpge->SetVisibility(true);
  vishpge->SetForceSolid(true);
  logicCryst->SetVisAttributes (vishpge);
  G4VisAttributes* visbewindow=new G4VisAttributes(G4Colour(0,1,0));
  visbewindow->SetVisibility(true);
  visbewindow->SetForceSolid(true);
  logicPatient->SetVisAttributes(visbewindow);
  G4VisAttributes* visdeadlayer=new G4VisAttributes(G4Colour(0,0,1));
  visdeadlayer->SetVisibility(true);
  visdeadlayer->SetForceSolid(true);
  logicDeadLayer->SetVisAttributes(visdeadlayer);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // declare crystal as a MultiFunctionalDetector scorer
  //
  G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
  SetSensitiveDetector("CrystalLV",cryst);


  // declare patient as a MultiFunctionalDetector scorer
  //
  G4MultiFunctionalDetector* patient = new G4MultiFunctionalDetector("patient");
  G4SDManager::GetSDMpointer()->AddNewDetector(patient);
  G4VPrimitiveScorer* primitiv2 = new G4PSDoseDeposit("dose");
  patient->RegisterPrimitive(primitiv2);
  SetSensitiveDetector("PatientLV",patient);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
