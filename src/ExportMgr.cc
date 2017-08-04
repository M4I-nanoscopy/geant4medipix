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
/// \file src/ExportMgr.cc
/// \brief Implementation of the ExportMgr class

// #include "G4GenericMessenger.hh"
#ifdef WITH_HDF5

#include "ExportMgr.hh"
#include "RunAction.hh"

#include "ExportHDF.hh"

#include "G4AutoLock.hh"
#include "locale.h"

ExportMgr *ExportMgr::instance = 0;

ExportMgr::ExportMgr()
{
  nbEvents = 0;
  filename = "";

  hdfExport = new ExportHDF(filename);
  //FIXME: is this still necessary?
  //exportMgr = new MPXExport();  // decide here which export to use

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExportMgr::~ExportMgr()
{
  delete hdfExport;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace
{
G4Mutex singletonMutex = G4MUTEX_INITIALIZER;
}

ExportMgr *ExportMgr::GetInstance()
{
  G4AutoLock l(&singletonMutex);
  if (!instance) {
      instance = new ExportMgr();
  }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace
{
G4Mutex AddDataMutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportMgr::AddData(DetectorHitsCollection *HitsCollection, G4int event, G4double energy)
{
  G4AutoLock l2(&AddDataMutex);
  // call export only when filename is set
  if (filename != "") {
    hdfExport->AddSingleEvents(HitsCollection, event);
    lastEvent = event;
    nbEvents++;
    if (nbEvents == 1000) {
      hdfExport->Write("/trajectories/", lastEvent, energy);
      nbEvents = 0;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ExportMgr::WriteData(G4double energy)
{
    if (filename != "") {
        hdfExport->Write("/trajectories/", lastEvent, energy);
        return lastEvent;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportMgr::SetHDFFilename(G4String name)
{
  ExportMgr *mgr = this->GetInstance();
  mgr->SetFilenameHDFexport(name);

  G4cout << "DEBUG: ExportMgrSetFilename: " << name << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace
{
G4Mutex SetFilenameMutex = G4MUTEX_INITIALIZER;
}

void ExportMgr::SetFilenameHDFexport(G4String name)
{
  G4AutoLock l2(&AddDataMutex);
  hdfExport->SetFilename(name);
  G4cout << "DEBUG: SetFilenameHDFexport: " << name << G4endl;

  filename = name;
}
#endif