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
/// \file src/EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"

#include "Run.hh"
#include "DetectorSD.hh"
#include "DetectorHit.hh"
#include "HistoManager.hh"

#include "DigitizerWeightField.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4GenericMessenger.hh"
#include "G4UnitsTable.hh"
#include "G4DigiManager.hh"
#include "DetectorConstructionBase.hh"


#include "Randomize.hh"
#include <iomanip>
#include <RunAction.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction() : G4UserEventAction(),
    fSensorHCID(-1),
    fPrintModulo(1000),
    fEnergyPerEvent(0)
{
    G4DigiManager *digiManager = G4DigiManager::GetDMpointer();
    DigitizerWeightField *wfDigitizer = new DigitizerWeightField("DigitizerWeightField");
    digiManager->AddNewModule(wfDigitizer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorHitsCollection *
EventAction::GetHitsCollection(G4int hcID,
                               const G4Event *event) const
{
    DetectorHitsCollection *hitsCollection
        = static_cast<DetectorHitsCollection *>(
              event->GetHCofThisEvent()->GetHC(hcID));

    if (! hitsCollection) {
        G4cerr << "Cannot access hitsCollection ID " << hcID << G4endl;
        exit(1);
    }

    return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::PrintEventStatistics(G4double sensorEdep,
                                       G4double sensorTrackLength) const
{
    // print event statistics
    G4cout
            << "   Sensor: total energy: "
            << std::setw(7) << G4BestUnit(sensorEdep, "Energy")
            << "       total track length: "
            << std::setw(7) << G4BestUnit(sensorTrackLength, "Length")
            << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *event)
{
    G4int eventID = event->GetEventID();
    if (eventID % fPrintModulo == 0 && eventID != 0) {
        G4cout << "\n---> Begin of event: " << eventID << G4endl;
    }
    // reset energy per event
    fEnergyPerEvent = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *event)
{
    // Print per event (modulo n)
    G4int eventID = event->GetEventID();
    if (eventID % fPrintModulo == 0 && eventID != 0) {
        G4cout << "---> End of event: " << eventID << G4endl;
    }

    G4RunManager *fRM = G4RunManager::GetRunManager();
    DetectorConstructionBase *myDet = (DetectorConstructionBase *)(fRM->GetUserDetectorConstruction());

    G4String digitizerName = myDet->GetDigitizerName();

    // call digitizer after every event
    G4DigiManager *digiManager = G4DigiManager::GetDMpointer();
    DigitizerWeightField *digiModule = (DigitizerWeightField *) (digiManager->FindDigitizerModule(digitizerName));
    digiModule->Digitize();

    if (fEnergyPerEvent > 0.) {
        //fill histogram with total energy per event
        G4AnalysisManager::Instance()->FillH1(2, fEnergyPerEvent);
    }
 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
