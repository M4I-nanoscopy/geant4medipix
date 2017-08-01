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
/// \file src/PrimaryGeneratorMessenger.hh
/// \brief Declaration of a generator messenger class

#include "PrimaryGeneratorMessenger.hh"

#include "G4RunManager.hh"

#include "G4GenericMessenger.hh"

#include "DetectorConstructionBase.hh"
#include "MpxDetector.hh"
#include "EventAction.hh"

#include "G4UIcommand.hh"
#include "G4SystemOfUnits.hh"

PrimaryGeneratorMessenger *PrimaryGeneratorMessenger::instance = 0;

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger() 
{
    // Default values
    nFrames = 1;
    step = 1;
    nPar = 1000;
    DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
    delete fMessenger;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace
{
G4Mutex singletonMutex = G4MUTEX_INITIALIZER;
}

PrimaryGeneratorMessenger *PrimaryGeneratorMessenger::GetInstance()
{
    G4AutoLock l(&singletonMutex);
    if (!instance) {
        instance = new PrimaryGeneratorMessenger();
    }
    return instance;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::DefineCommands()
{
    // create a folder and command for dacscans 
    fMessenger = new G4GenericMessenger(this, "/Medipix/gun/", "Control the run");
    fMessenger->DeclareMethod("beamOn", &PrimaryGeneratorMessenger::beamOn, "Beam On.");
    fMessenger->DeclareMethod("frames", &PrimaryGeneratorMessenger::setNumberOfFrames, "Set the number of frames");
    // define the dacscan folder
    fMessenger->SetDirectory("/Medipix/dacscan/");
    fMessenger->DeclareMethod("particles", &PrimaryGeneratorMessenger::SetNumberOfParticles, "Set the number of particles in one step");
    fMessenger->DeclareMethod("step", &PrimaryGeneratorMessenger::SetStepSize, "Set the step size of the dacscan");
    fMessenger->DeclareMethod("run", &PrimaryGeneratorMessenger::Dacscan, "run a dacscan [start] [stop]");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorMessenger::beamOn(G4int np)
{

    G4cout << "Total number of frames: " << nFrames << G4endl;

    //Get RunManager, DetectorConstruction and Detector
    G4RunManager *runManager = G4RunManager::GetRunManager();
    DetectorConstructionBase *myDet = (DetectorConstructionBase *)(runManager->GetUserDetectorConstruction());
    MpxDetector *mpxDet = MpxDetector::GetInstance();

    G4int detectorType;
    detectorType = myDet->GetDetectorType();
    G4cout << "Detector type " << detectorType << G4endl;
    // 0 = Medipix3RX
    // 1 = Timepix
    // 2 = Dosepix

    G4String base = myDet->GetOutputFilename();
    G4String thStr;

    for (int i = 0; i < nFrames; i++) {
        //reset pixel matrix
        mpxDet->ResetPixelMatrix();
	mpxDet->UpdateSettings();
        // set threshold for medipix and Timepix detector
        if (detectorType == 0) {
            G4double th1, th2, th3, th4;
            th1 = mpxDet->GetMpxThreshold1();
            th2 = mpxDet->GetMpxThreshold2();
            th3 = mpxDet->GetMpxThreshold3();
            th4 = mpxDet->GetMpxThreshold4();
            G4cout << "Frame: " << i << " Getting ready to fire " << np << " particles. Thresholds are at "
                   << th1 << " " << th2 << " " << th3 << " and " << th4 << " keV" << G4endl;
                   
        } else if (detectorType == 1) {
            G4double th;
            th = mpxDet->GetTpxThreshold();
            G4cout << "Frame: " << i << " Getting ready to fire " << np << " particles. Threshold is at " << th << " keV" << G4endl;
        }
        thStr = G4UIcommand::ConvertToString(i);
	//if (base != "")
	//    myDet->SetOutputFilename(base + "_" + thStr);
        runManager->BeamOn(np);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorMessenger::Dacscan(G4double start, G4double stop)
{

    G4RunManager *runManager = G4RunManager::GetRunManager();
    DetectorConstructionBase *myDet = (DetectorConstructionBase *)(runManager->GetUserDetectorConstruction());
    MpxDetector *mpxDet = MpxDetector::GetInstance();
    G4String base = myDet->GetOutputFilename();
    G4String thStr;

    G4double th = start;
    // filename for dacsteps
    while (th < stop) {
        thStr = G4UIcommand::ConvertToString(th * 1000);
        mpxDet->SetTpxThreshold(th * keV);
        myDet->SetOutputFilename(base + "_Threshold_" + thStr + "eV.mpx");
        this->beamOn(nPar);
        th += step;
    }
    myDet->SetOutputFilename(base);

}
