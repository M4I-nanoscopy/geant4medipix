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
//
/// \file src/Run.cc
/// \brief Implementation of the Run class

#include "Run.hh"

#include "DetectorHit.hh"
#include "G4SDManager.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run() : G4Run(),
fEdeposit(0.),
fEdeposit2(0.)
{
#ifdef WITH_HDF5
    mgr = ExportMgr::GetInstance();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent(const G4Event *event)
{
#ifdef WITH_HDF5
    DetectorHitsCollection *HitsCollection = new DetectorHitsCollection;
    HitsCollection = static_cast<DetectorHitsCollection *>(event->GetHCofThisEvent()->GetHC(0));

    lastEvent = event->GetEventID();

    if (HitsCollection->GetSize() != 0) {
        mgr->AddData(HitsCollection, lastEvent);
    }
#endif
}

void Run::Merge(const G4Run* run)
{
  const Run* localRun =  static_cast<const Run*>(run);
//   fEdeposit  += localRun->fEdeposit;
//   fEdeposit2 += localRun->fEdeposit2;
  
  G4Run::Merge(run);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

