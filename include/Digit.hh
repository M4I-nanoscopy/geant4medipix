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
/// \file Digit.hh
/// \brief Definition of the MpxDigit class

#ifndef MpxDigit_h 
#define MpxDigit_h 1 

#include "G4VDigi.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"

class Digit : public G4VDigi
{
public:
  Digit();
  ~Digit();  
 
  inline G4int  GetColumn(){ return column;}
  inline void   SetColumn(G4int col){ column = col;}
  
  inline G4int  GetLine(){ return line;}
  inline void   SetLine(G4int lin){ line = lin;}
  
  inline G4int  GetEvent(){ return event;}
  inline void   SetEvent(G4int ev){event = ev;}
  
  inline G4double   GetEnergy(){ return energy;}
  inline void       SetEnergy(G4double en){ energy = en;}

  inline G4double   GetToT(){ return tot;}
  inline void       SetToT(G4double ToT){ tot = ToT;}

  inline G4double   GetToA(){ return toa;}
  inline void       SetToA(G4double ToA){ toa = ToA;}
  
  inline G4double*  GetPreAmpResp(){ return preAmpResp;}
  inline void       SetPreAmpResp(G4double* preAmp){ preAmpResp = preAmp;}
  
private:
    G4double        energy;     //collected energy
    G4double        tot;        //ToT value
    G4double        toa;        //ToA value
    G4int           column;     //column of Pixel
    G4int           line;       //line of Pixel
    G4int           event;      //event id for spare particle based format
    //FIXME needed?
    G4double*       preAmpResp = NULL; //response curve
};

typedef G4TDigiCollection<Digit> MpxDigitCollection;

#endif
