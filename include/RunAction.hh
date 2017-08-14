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
/// \file include/RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "MpxDetector.hh"
#include "globals.hh"
#include "H5Cpp.h"
#include "ExportMgr.hh"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
class G4Run;
class G4Timer;
class HistoManager;
class G4GenericMessenger;

// Run action class


class RunAction : public G4UserRunAction

{
public:
    RunAction();
    virtual ~RunAction();

    virtual G4Run *GenerateRun();
    virtual void   BeginOfRunAction(const G4Run *);
    virtual void   EndOfRunAction(const G4Run *);
    ExportMgr *getExportManager() const;
private:
    G4Timer *timer;
    MpxDetector* detector;
    HistoManager* histoManager;
    G4GenericMessenger* fMessenger;
    ExportMgr *exportManager;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

