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

#include <ExportRaw.hh>
#include <ExportHDF.hh>
#include "ExportMgr.hh"

ExportMgr::ExportMgr()
{
  nbEvents = 0;

  rawExport = new ExportHDF();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExportMgr::~ExportMgr()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace
{
G4Mutex AddDataMutex = G4MUTEX_INITIALIZER;
}

bool writeAttributes = false;
const int PIXELS_CHUNK_SIZE = 100;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportMgr::AddData(DetectorHitsCollection *HitsCollection, MpxDigitCollection *DigitCollection) {
    G4AutoLock l2(&AddDataMutex);

    // call export only when filename is set
    if (!filename.empty()) {
        rawExport->AddSingleEvents(HitsCollection);
        rawExport->AddSingleDigits(DigitCollection);

        nbEvents++;
        // We used to write pixels and trajectories in chunks to the H5 file. This turned out to be risky in terms
        // of segfaults and file locks. There is no real drawback to do all writing at the end (memory is big enough to
        // hold all results)
        if (nbEvents == (1)) {
            //rawExport->Write();
            //rawExport->WritePixels();
            nbEvents = 0;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportMgr::WriteData()
{
    G4AutoLock l2(&AddDataMutex);

    if (!filename.empty()) {
        //rawExport->Write();
        //rawExport->WritePixels();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportMgr::CreateDataFile() {
    G4AutoLock l2(&AddDataMutex);

    if (!filename.empty()) {
        rawExport->CreateOutputFile();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportMgr::SetHDFFilename(G4String name)
{
    rawExport->SetFilename(name);
    filename = name;
}

void ExportMgr::SetAttributes() {
    G4AutoLock l2(&AddDataMutex);

    // Only write attributes once by using global writeAttributes bool
    if (filename != "" && !writeAttributes) {
        rawExport->SetAttributes();
        writeAttributes = true;
    }
}


#endif