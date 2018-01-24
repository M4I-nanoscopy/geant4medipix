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
/// \file src/DetectorSD.cc
/// \brief Implementation of the DetectorSD class

#include "DetectorSD.hh"

#include "DetectorConstructionBase.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorSD::DetectorSD(
    const G4String &name,
    const G4String &hitsCollectionName,
    G4int nofPixels)
    : G4VSensitiveDetector(name),
      fHitsCollection(0),
      fNofCells(1)
{
    collectionName.insert(hitsCollectionName);
    nPixels 	= nofPixels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorSD::~DetectorSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorSD::Initialize(G4HCofThisEvent *hce)
{
    //G4int hcID = GetCollectionID(0);

    // Create hits collection
    fHitsCollection
        = new DetectorHitsCollection(SensitiveDetectorName, collectionName[0]);

    // Add this collection in hce
    G4int hcID
        = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(hcID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Alpha.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"

G4bool DetectorSD::ProcessHits(G4Step *step, G4TouchableHistory *)
{
    // energy deposit
    G4double edep = step->GetTotalEnergyDeposit();
    G4ThreeVector pos = step->GetTrack()->GetPosition();

    // get particle and create ID
    G4String particle = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    G4int particleID = 0;

    if (particle == "gamma") {
        particleID = 0;
    } else if (particle == "e-") {
        particleID = 1;
    } else if (particle == "alpha") {
        particleID = 3;
    } else if (particle == "proton") { // anti_proton?
        particleID = 2;
    } else if (particle == "e+") {
        particleID = -1;
    } else {
        //particle ID not defined
        particleID = 9999;
        //neutron,anti_proton
    }

    // track global time
    G4double time = 0;
    time = step->GetPostStepPoint()->GetLocalTime();

    // step length
    G4double stepLength = 0.;
    // get step length after energy deposit
    if (step->GetTrack()->GetDefinition()->GetPDGCharge() != 0.) {
        stepLength = step->GetStepLength();
    } 
    // add hit if energy deposit > 0 and step particle came to rest
    if (edep == 0. && stepLength == 0.) {
        return false;
    } else {
        G4TouchableHistory *touchable
            = (G4TouchableHistory *)(step->GetPreStepPoint()->GetTouchable());

        // Get Detector cell id
        G4int lineNumber = touchable->GetReplicaNumber(1);
        G4int columnNumber = touchable->GetReplicaNumber(0);

        //get currentEvent
        G4RunManager *fRM = G4RunManager::GetRunManager();
        const G4Event *currentEvent = fRM->GetCurrentEvent();
        eventNb = currentEvent->GetEventID();

        DetectorHit *hit = new DetectorHit();
        
        // add hit information to hits collection
        hit->Add(edep, stepLength, pos, columnNumber, lineNumber, eventNb, particleID, time);
        fHitsCollection->insert(hit);

        return true;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorSD::EndOfEvent(G4HCofThisEvent *)
{
    if (verboseLevel > 1) {
        G4int nofHits = fHitsCollection->entries();
        G4cout << "\n-------->Hits Collection: in this event there are " << nofHits
               << " hits in the detector: " << G4endl;
        for (G4int i = 0; i < nofHits; i++)(*fHitsCollection)[i]->Print();
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
