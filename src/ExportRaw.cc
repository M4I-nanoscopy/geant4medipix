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
// $Id: ExportRaw.cc 71060 2013-06-10 15:03:19Z gcosmo $
//
/// \file src/ExportRaw.cc
/// \brief Implementation of the ExportRaw class


#include "ExportRaw.hh"

#include "DetectorConstructionBase.hh"

#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"

#include <G4GenericMessenger.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExportRaw::ExportRaw()
{
    filename    = "g4medipix.raw";
}

namespace
{
    G4Mutex WriteMutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExportRaw::AddSingleEvents(DetectorHitsCollection *HitsCollection)
{
    /*for (G4int i = 0; i < (G4int) HitsCollection->GetSize(); i++) {
        // Hard copy of object
        //if (i > 0 && (*HitsCollection)[i]->GetTime() < (*HitsCollection)[i-1]->GetTime()) {
        //    return;
        //}
        DetectorHit *hitDetector = (*HitsCollection) [i];
        auto *hitCopy = new DetectorHit(*hitDetector);

        HitsCollectionCopy->insert(hitCopy);
    }*/
    ExportRaw::Write(HitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExportRaw::AddSingleDigits(MpxDigitCollection *DigitCollection)
{
    /*for (G4int i = 0; i < (G4int) DigitCollection->GetSize(); i++) {
        // Hard copy of object
        Digit *digit = (*DigitCollection) [i];
        DigitCollectionCopy->insert(new Digit(*digit));
    }*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportRaw::Write(DetectorHitsCollection *hc) {
    size_t LENGTH = hc->GetSize();
    size_t MAX_TRAJ = 1000;

    G4cout << "Writing trajectories output per event to raw. Number of hits: " << LENGTH << G4endl;

    if ( LENGTH == 0) {
        return;
    }

    // structure  and dataset
    struct s1_t {
        G4double x;
        G4double y;
        G4double z;
        G4double energy;
    };

    // Handles


    auto *det = (DetectorConstructionBase *)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    if (! det->storeTraj) {
        return;
    }

    // Mutex lock
    G4AutoLock autoLock(&WriteMutex);

    // Get file
    int file = GetOutputFile();

    // Array of structs to store
    auto *s1 = new s1_t[MAX_TRAJ];
    memset(s1, 0, MAX_TRAJ * sizeof(struct s1_t));

    G4int ev = (*hc)[0]->GetEvent();
    G4int temp = 0;

    for (size_t i = 0; i < LENGTH; i++) {

        DetectorHit *sensorHit = (*hc)[i];

        if (sensorHit->GetEvent() != ev || i == LENGTH - 1) {
            if (i == LENGTH - 1) {
                s1[temp].x = sensorHit->GetPosition().x() / nm;
                s1[temp].y = sensorHit->GetPosition().y() / nm;
                s1[temp].z = sensorHit->GetPosition().z() / nm;
                s1[temp].energy = sensorHit->GetEdep() / keV;
            }

            memset(s1, 0, MAX_TRAJ * sizeof(struct s1_t));
            temp = 0;

            ev = sensorHit->GetEvent();

        }

        s1[temp].x = sensorHit->GetPosition().x() / nm;
        s1[temp].y = sensorHit->GetPosition().y() / nm;
        s1[temp].z = sensorHit->GetPosition().z() / nm;
        s1[temp].energy = sensorHit->GetEdep() / keV;
        temp++;
    }

    // Close output (first close output, then destroy memory)
    CloseOutputFile();

    // New HitsCollection
    delete[] s1;
}

void ExportRaw::WritePixels(MpxDigitCollection *dc) {
    size_t number_digits = dc->GetSize();

    if (number_digits == 0) {
        return;
    }

    auto *det = (DetectorConstructionBase *)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    auto nb = (size_t) det->GetNbPixels();

    // Reserve space
    auto *pixels = new G4double[2 * nb * nb]{0};

    // Mutex lock
    G4AutoLock autoLock(&WriteMutex);

    // Handles
    int file = GetOutputFile();

    // Get first event
    G4int event = (*dc)[0]->GetEvent();

    // Iterate over all entries in the list
    for (size_t i = 0; i < dc->GetSize(); i++) {

        Digit *d = (*dc)[i];

        if (d->GetEvent() == 0) {
            G4cout << "Found a digit with event 0. Throwing it away, not trusting it." << G4endl;
            continue;
        }

        // Digit collection start at 1 (hence -1)
        size_t x = (d->GetColumn() - 1)*nb + (d->GetLine() - 1);
        size_t y = nb*nb + (d->GetColumn() - 1)*nb + (d->GetLine() - 1);

        if ( event != d->GetEvent() || i == dc->GetSize() - 1) {

            if ( i == dc->GetSize() - 1 ) {
                pixels[x] = d->GetToT();
                pixels[y] = d->GetToA();
            }

            G4cout << "Writing sparse pixels output per event " << event << " to raw. Digits: " << i << G4endl;

            // Open dataset

            // Set all ToA values to 0 and all ToT values to 0
            memset(pixels, 0, 2 * nb * nb * sizeof(G4double));

            event = d->GetEvent();
        }

        pixels[x] = d->GetToT();
        pixels[y] = d->GetToA();
    }

    // Close output (first close output, then destroy memory)
    CloseOutputFile();

    // Start new DigitCollection
    delete[] pixels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportRaw::SetFilename(G4String name)
{
// #ifdef G4MULTITHREADED
    // std::ostringstream os;
    // os << G4Threading::G4GetThreadId();
    //name.append("_t");
    //name.append(os.str());
// #endif

    filename = name;
}

void ExportRaw::SetAttributes() {
    auto *det = (DetectorConstructionBase *)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    G4double height = det->GetSensorThickness() / nm;
    G4String mat = det->GetSensorMaterial()->GetName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportRaw::CreateOutputFile() {
    G4AutoLock autoLock(&WriteMutex);

    G4cout << "Creating raw output file " << filename.c_str() << G4endl;

    // Create file and dataset

    // ExportRaw::SetAttributes();

    G4cout << "Closed" << G4endl;
}

int ExportRaw::GetOutputFile() {
    int file = 0;

    return file;
}

void ExportRaw::CloseOutputFile() {

    G4cout << "Closed" << G4endl;
}

