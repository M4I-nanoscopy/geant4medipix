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
/// \file src/Digit.cc
/// \brief Implementation of the MpxDigit class

#include "Digit.hh"

G4ThreadLocal G4Allocator<Digit> *DigitAllocator = 0;

Digit::Digit() :
    energy(0),
    tot(0),
    toa(0),
    column(0),
    line(0),
    event(0)
{
}

Digit::~Digit()
{}

Digit::Digit(const Digit& right)
        :G4VDigi()
{
    energy = right.energy;
    toa = right.toa;
    tot = right.tot;
    column = right.column;
    line = right.line;
    event = right.event;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const Digit& Digit::operator=(const Digit& right)
{
    energy = right.energy;
    toa = right.toa;
    tot = right.tot;
    column = right.column;
    line = right.line;
    event = right.event;
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int Digit::operator==(const Digit& right) const
{
    return (
            (energy == right.energy) &&
            (toa == right.toa) &&
            (tot == right.tot) &&
            (column == right.column) &&
            (line == right.line) &&
            (event == right.event)
    );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Digit::Draw()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Digit::Print()
{;}