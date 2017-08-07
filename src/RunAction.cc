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
/// \file src/RunAction.cc
/// \brief Implementation of the RunAction class

#include <DetectorConstructionBase.hh>
#include <G4DigiManager.hh>
#include <DigitizerWeightField.hh>
#include "RunAction.hh"
#include "DetectorHit.hh"
#include "DetectorSD.hh"
#include "HistoManager.hh"
#include "Run.hh"

#ifdef WITH_HDF5
#include "ExportMgr.hh"
#endif
#include "G4Timer.hh"


#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"

//static G4ThreadLocal Messenger *fmessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
    : G4UserRunAction(),
      histoManager(0)
{
#ifdef WITH_HDF5
    G4String name;
    fMessenger = new G4GenericMessenger(this, "/Output/", "Define all file names here");
    fMessenger->DeclareMethod("edep", &ExportMgr::SetHDFFilename, "The filename for the HDF5 export.").SetParameterName(name, true).SetStates(G4State_Idle, G4State_PreInit);
#endif
    histoManager = new HistoManager();
    timer = new G4Timer();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
#ifdef WITH_HDF5
    delete fMessenger; // TODO check if ok to kill here
#endif
    delete histoManager;

}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run *RunAction::GenerateRun()
{
    return new Run();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run *)
{
    if (isMaster)
        timer->Start();

    // Build the digitizer
    G4RunManager *fRM = G4RunManager::GetRunManager();
    DetectorConstructionBase *myDet = (DetectorConstructionBase *)(fRM->GetUserDetectorConstruction());
    G4String digitizerName = myDet->GetDigitizerName();
    G4DigiManager *digiManager = G4DigiManager::GetDMpointer();

    if (digitizerName == "WFDigitizer") {
        DigitizerWeightField *wfDigitizer = new DigitizerWeightField("WFDigitizer");
        digiManager->AddNewModule(wfDigitizer);
    }

    // create histogram file
    G4AnalysisManager *analysisManager =  G4AnalysisManager::Instance();
    if (analysisManager->IsActive()) {
        analysisManager->OpenFile();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run *aRun)
{
  if (isMaster) {
      timer->Stop();

      G4cout << "Run " << aRun->GetRunID() << " with "
              << aRun->GetNumberOfEventToBeProcessed()
              << " particles processed in " << timer->GetRealElapsed()
              << "s" << G4endl << G4endl;

      detector = MpxDetector::GetInstance();
      detector->WriteSparse();
      detector->WriteFrame();
      // detector->WriteSimulationSettings();

#ifdef WITH_HDF5
        //write Data to HDF5 file and delet Export Manager
        ExportMgr *mgr = ExportMgr::GetInstance();
      mgr->WriteData();
#endif
    }
    // write histogram files
    // Analysis manager takes care of threads and joins files!
    G4AnalysisManager *analysisManager =  G4AnalysisManager::Instance();
    if (analysisManager->IsActive()) {
        analysisManager->Write();
        analysisManager->CloseFile();
    }
}




