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

  hdfExport = new ExportHDF();
  //FIXME: is this still necessary?
  //exportMgr = new MPXExport();  // decide here which export to use

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExportMgr::~ExportMgr()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExportMgr *ExportMgr::GetInstance()
{
    G4bool isMaster = ! G4Threading::IsWorkerThread();

    if (!instance && isMaster) {
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

void ExportMgr::AddData(DetectorHitsCollection *HitsCollection, MpxDigitCollection *DigitCollection, G4int event) {
    G4AutoLock l2(&AddDataMutex);

    // call export only when filename is set
    if (filename != "") {
        hdfExport->AddSingleEvents(HitsCollection);
        hdfExport->AddSingleDigits(DigitCollection);

        lastEvent = event;
        nbEvents++;
        if (nbEvents == PIXELS_CHUNK_SIZE) {
            hdfExport->Write("/trajectories/", lastEvent);
            hdfExport->WritePixels();
            nbEvents = 0;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportMgr::WriteData()
{
    G4AutoLock l2(&AddDataMutex);

    if (filename != "") {
        hdfExport->Write("/trajectories/", lastEvent);
        hdfExport->WritePixels();
    }
}

void ExportMgr::CreateDataFile() {
    G4AutoLock l2(&AddDataMutex);

    hdfExport->CreateOutputFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportMgr::SetHDFFilename(G4String name)
{
    G4cout << "DEBUG: SetFilenameHDFexport: " << name << G4endl;
    hdfExport->SetFilename(name);
    filename = name;
}

#endif