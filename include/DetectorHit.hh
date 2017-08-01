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
/// \file DetectorHit.hh
/// \brief Definition of the DetectorHit class

#ifndef DetectorHit_h
#define DetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

// Detectorimeter hit class
//
// It defines data members to store the the energy deposit and track lengths
// of charged particles in a selected volume:
// - fEdep, fTrackLength

class DetectorHit : public G4VHit
{
public:
    DetectorHit();
    DetectorHit(const DetectorHit &);
    virtual ~DetectorHit();

    // operators
    const DetectorHit &operator=(const DetectorHit &);
    G4int operator==(const DetectorHit &) const;

    inline void *operator new(size_t);
    inline void  operator delete(void *);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void Add(G4double de, G4double dl, G4ThreeVector pos, G4int column, G4int line, G4int event, G4int particleID, G4double time);
    void SetEdep(G4double);
    
    // get methods
    G4double GetEdep() const;
    G4double GetTrackLength() const;
    G4ThreeVector GetPosition() const;
    G4int GetColumn() const;
    G4int GetLine() const;
    G4int GetEvent() const;
    G4int GetParticleID() const;
    G4double GetTime() const;

private:
    G4double        fEdep;        	//< Energy deposit in the sensitive volume
    G4double        fTrackLength; 	//< Track length in the  sensitive volume
    G4ThreeVector   fPosition;	//
    G4int           fColumn;      	//column of Pixel
    G4int           fLine;	  	//line of Pixel
    G4int           fEvent;		//eventNb
    G4int           fParticleID;
    G4double        fTime;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<DetectorHit> DetectorHitsCollection;

extern G4ThreadLocal G4Allocator<DetectorHit> *DetectorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void *DetectorHit::operator new(size_t)
{
    if (!DetectorHitAllocator)
        DetectorHitAllocator = new G4Allocator<DetectorHit>;
    void *hit;
    hit = (void *) DetectorHitAllocator->MallocSingle();
    return hit;
}

inline void DetectorHit::operator delete(void *hit)
{
    if (!DetectorHitAllocator)
        DetectorHitAllocator = new G4Allocator<DetectorHit>;
    DetectorHitAllocator->FreeSingle((DetectorHit *) hit);
}

inline void DetectorHit::Add(G4double de, G4double dl, G4ThreeVector pos,  G4int column, G4int line, G4int event, G4int particleID, G4double time)
{
    fEdep 	+= de;
    fTrackLength += dl;
    fPosition 	+= pos;
    fColumn 	= column;
    fLine 	= line;
    fEvent 	= event;
    fParticleID = particleID;
    fTime = time;
}

inline G4int DetectorHit::GetColumn() const
{
    return fColumn;
}

inline G4int DetectorHit::GetLine() const
{
    return fLine;
}

inline G4double DetectorHit::GetEdep() const
{
    return fEdep;
}

inline void DetectorHit::SetEdep(G4double de)
{
    fEdep = de;
}

inline G4double DetectorHit::GetTrackLength() const
{
    return fTrackLength;
}

inline G4ThreeVector DetectorHit::GetPosition() const
{
    return fPosition;
}
inline G4int DetectorHit::GetEvent() const
{
    return fEvent;
}
inline G4int DetectorHit::GetParticleID() const
{
    return fParticleID;
}
inline G4double DetectorHit::GetTime() const
{
    return fTime;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
