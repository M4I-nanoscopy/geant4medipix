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
/// \file src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstructionBase.hh"

#include "DetectorMessenger.hh"
#include "DetectorSD.hh"

#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4AffineTransform.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4GenericMessenger.hh"
#include "G4AssemblyVolume.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GenericMessenger *DetectorConstructionBase::fMessenger = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstructionBase::DetectorConstructionBase() : G4VUserDetectorConstruction(),
    fCheckOverlaps(true),
    detectorType(0),
    csmMode(false),
    filter(false),
    collimator(false)
{
    // some default values for now
    nPixel = 10;
    pixelSize = 110 * um;
    sensorThickness = 300 * um;
    filterThickness = 3 * um;
    filterZ = .5 * mm;

    //Bumb bonds
    bumps = false;
    bumpRadii = 50 * um;
    bumpHeight = 50 * um;

    //electronics chip
    chip = false;
    pcb = false;

    //Output
    optFname = "";
    sparseOptFname = "";
    configFilename = "";

    //Collimator
    collimatorThickness = 500 * um;
    colRadii = pixelSize / 2. / 2.;

    block = false;
    actX = 14080 * um;
    actY = 14080 * um;
    actZ = 500 * um;
    actPosZ = -1. * mm;
    
    //rotation of sensor and other solids in world volume 
    fRotation.rotateX(0.*deg);
    fRotation.rotateY(0.*deg);
    fRotation.rotateZ(0.*deg);

    CalculateGeometry();

    // Define materials
    DefineMaterials();
    worldMaterial = G4Material::GetMaterial("G4_Galactic");
    sensorMaterial = G4Material::GetMaterial("G4_Si");
    filterMaterial = G4Material::GetMaterial("G4_Al");
    collimatorMaterial = G4Material::GetMaterial("G4_Fe");
    bumpMaterial = G4Material::GetMaterial("G4_Cu");
    // set for Medipix3RX in case gain is empty
    feedbackCapacitance = 8e-15;
    
    fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstructionBase::~DetectorConstructionBase()
{
    delete fMessenger;
    delete fDetectorMessenger;

    fMessenger = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstructionBase::CalculateGeometry()
{
    sensorSizeXY = nPixel * pixelSize;

    G4cout << "NoPixels: " << nPixel << G4endl;
    G4cout << "sensorThickness: " << sensorThickness << G4endl;
    G4cout << "sensorSizeXY: " << sensorSizeXY << G4endl;

//    G4cout << "sensorMaterial: " << sensorMaterial << G4endl;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstructionBase::Construct()
{
    //clean out old geometry before updating
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    // Define volumes
    return DefineVolumes();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstructionBase::DefineMaterials()
{
    // standard material defined using NIST Manager
    G4NistManager *nistManager = G4NistManager::Instance();

    nistManager->FindOrBuildMaterial("G4_Galactic");
    nistManager->FindOrBuildMaterial("G4_Si");
    nistManager->FindOrBuildMaterial("G4_CADMIUM_TELLURIDE");
    nistManager->FindOrBuildMaterial("G4_Ni");
    nistManager->FindOrBuildMaterial("G4_Ti");
    nistManager->FindOrBuildMaterial("G4_Pr");
    nistManager->FindOrBuildMaterial("G4_Ag");
    nistManager->FindOrBuildMaterial("G4_Au");
    nistManager->FindOrBuildMaterial("G4_Cu");
    nistManager->FindOrBuildMaterial("G4_Pd");
    nistManager->FindOrBuildMaterial("G4_Zr");
    nistManager->FindOrBuildMaterial("G4_Al");
    nistManager->FindOrBuildMaterial("G4_Fe");
    nistManager->FindOrBuildMaterial("G4_In");
    nistManager->FindOrBuildMaterial("G4_Pb");
    nistManager->FindOrBuildMaterial("G4_W");
    nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    nistManager->FindOrBuildMaterial("G4_AIR");
    nistManager->FindOrBuildMaterial("G4_Np");
    nistManager->FindOrBuildMaterial("G4_GALLIUM_ARSENIDE");

    //Additional materials
    G4double z, a, fractionmass, density;
    G4String name, symbol;
    G4int ncomponents;

    a = 118.71 * g / mole;
    G4Element *elSn  = new G4Element(name = "Tin", symbol = "Sn" , z = 50., a);

    a = 207.2 * g / mole;
    G4Element *elPb  = new G4Element(name = "Lead"  , symbol = "Pb" , z = 82., a);

    density = (1 / (0.6 / 7.365 + 0.4 / 11.34)) * g / cm3;
    G4Material *SnPb = new G4Material(name = "Tin/Lead", density, ncomponents = 2);
    SnPb->AddElement(elSn, fractionmass = 50 * perCent);
    SnPb->AddElement(elPb, fractionmass = 50 * perCent);

    // Print available materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstructionBase::ConstructSDandField()
{
    // Sensitive detectors
    DetectorSD *sensorSD;

    // check if already created
    sensorSD = dynamic_cast<DetectorSD *>(G4SDManager::GetSDMpointer()->FindSensitiveDetector("SensorSD", false));

    if (!sensorSD) {
        // create a new sensitive detector
        sensorSD = new DetectorSD("SensorSD", "SensorHitsCollection", nPixel * nPixel);
    }

    SetSensitiveDetector("pixel_cell", sensorSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetWorldMaterial(const G4String &material)
{
  G4Material *pMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pMaterial)
    worldMaterial = pMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetSensorThickness(G4double thickness)
{
    if (thickness <= DBL_MIN) {
        G4cout << "\n ---> warming from SetDetectorThickness: thicknes "
               << thickness << " out of range." << G4endl;
        return;
    }
    sensorThickness = thickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetSensorMaterial(const G4String &material)
{
    G4Material *pMaterial =
        G4NistManager::Instance()->FindOrBuildMaterial(material);
    if (pMaterial)
        sensorMaterial = pMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetFilterMaterial(const G4String &material)
{
    G4Material *pMaterial =
        G4NistManager::Instance()->FindOrBuildMaterial(material);
    if (pMaterial)
        filterMaterial = pMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetNbPixels(G4int nbpixels)
{
    if (nbpixels < 1) {
        G4cout << "\n ---> warning from SetNumberPixels: "
               << nbpixels << " out of range." << G4endl;
        return;
    }
    nPixel = nbpixels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetPixelSize(G4double pixsize)
{

    G4RunManager::GetRunManager()->GeometryHasBeenModified();

    if (pixsize <= 0.0) {
        G4cout << "\n ---> warning from SetPixelSize: "
               << pixsize << " out of range." << G4endl;
        return;
    }
    pixelSize = pixsize;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstructionBase::SetDetectorType(G4String dname)
{
    if (dname == "Medipix3RX") {
        detectorType = 0;
    } else if (dname == "Timepix") {
        detectorType = 1;
    } else if (dname == "Dosepix") {
        detectorType = 2;
    } else if (dname == "Timepix3") {
        detectorType = 3;
    }
}

void DetectorConstructionBase::SetDetectorGain(G4String dname)
{
    if (dname == "SLGM") {
        feedbackCapacitance = 28e-15;
    } else if (dname == "LGM") {
        feedbackCapacitance = 21e-15;
    } else if (dname == "HGM") {
        feedbackCapacitance = 14e-15;
    } else if (dname == "SHGM") {
        feedbackCapacitance = 7e-15;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetTpxMode(G4String mode)
{
//FIXME: do we still need this function since we write everything at once?
  if (mode == "PC") {
    tpxMode = 0;
  } else if (mode == "ToT") {
    tpxMode = 1;
  } else if (mode == "ToA") {
    tpxMode = 2;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetRotation(G4double xangle)
{
  if (xangle < 90 || xangle > 90){
    fRotation.rotateX(xangle*deg);
    fRotation.rotateY(0*deg);
    fRotation.rotateZ(0*deg);
  }
  else {
    G4cout << "\n ---> warning from SetRotation: angle "
    << xangle << " out of range." << G4endl;
    return;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetFilterThickness(G4double thickness)
{
    if (thickness <= DBL_MIN) {
        G4cout << "\n ---> warning from SetDetectorThickness: thickness "
               << thickness << " out of range." << G4endl;
        return;
    }
    filterThickness = thickness;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetFilterZ(G4double z)
{
    if (z <= DBL_MIN) {
        G4cout << "\n dgb78 ---> warning from SetFilterZ: Z "
               << z << " out of range." << G4endl;
        return;
    }
    G4cout << "dbg77: " << z << G4endl;
    filterZ = z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetFilterRotation(G4double xangle){
    if (xangle < 90 || xangle > 90){
      filterRotation.rotateX(xangle*deg);
      filterRotation.rotateY(0*deg);
      filterRotation.rotateZ(0*deg);
    }
    else {
      G4cout << "\n ---> warning from SetFilterRotation: angle "
      << xangle << " out of range." << G4endl;
      return;
    }
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// bump configuration as inline functions: see header

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetBumpMaterial(const G4String &material)
{
    G4Material *pMaterial =
        G4NistManager::Instance()->FindOrBuildMaterial(material);
    if (pMaterial)
        bumpMaterial = pMaterial;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// chip, output filename, sparse as inline functions: see header

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::SetCollimatorMaterial(const G4String &material)
{
    G4Material *pMaterial =
        G4NistManager::Instance()->FindOrBuildMaterial(material);
    if (pMaterial)
        collimatorMaterial = pMaterial;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// collimator and block as inline functions: see header

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstructionBase::UpdateGeometry()
{
    //reset Geometry
    G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
    //reset SensitiveDetector
    ConstructSDandField();
}
