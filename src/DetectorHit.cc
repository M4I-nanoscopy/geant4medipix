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
/// \file src/DetectorHit.cc
/// \brief Implementation of the DetectorHit class

#include "DetectorHit.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<DetectorHit> *DetectorHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorHit::DetectorHit() : G4VHit(),
      fEdep(0.),
      fTrackLength(0.),
      fPosition(0., 0., 0.),
      fColumn(0),
      fLine(0),
      fEvent(0),
      fParticleID(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorHit::~DetectorHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorHit::DetectorHit(const DetectorHit &right)
    : G4VHit()
{
    fEdep        = right.fEdep;
    fTrackLength = right.fTrackLength;
    fPosition    = right.fPosition;
    fColumn      = right.fColumn;
    fLine        = right.fLine;
    fEvent  	 = right.fEvent;
    fParticleID  = right.fParticleID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const DetectorHit &DetectorHit::operator=(const DetectorHit &right)
{
    fEdep        = right.fEdep;
    fTrackLength = right.fTrackLength;
    fPosition    = right.fPosition;
    fColumn      = right.fColumn;
    fLine   	 = right.fLine;
    fEvent  	 = right.fEvent;
    fParticleID  = right.fParticleID;

    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DetectorHit::operator==(const DetectorHit &right) const
{
    return (this == &right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorHit::Print()
{
    G4cout
            << "Edep: "
            << std::setw(7) << G4BestUnit(fEdep, "Energy")
            << " track length: "
            << std::setw(7) << G4BestUnit(fTrackLength, "Length")
            << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
