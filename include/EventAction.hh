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
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"

#include "DetectorHit.hh"
#include "globals.hh"

class G4GenericMessenger;


class G4Track;
// Event action class
class EventAction : public G4UserEventAction
{
public:
    EventAction();
    virtual ~EventAction();

    virtual void  BeginOfEventAction(const G4Event *);
    virtual void  EndOfEventAction(const G4Event *);

    /**
     * Print every X event
     * \param value any positive integer
     */
    inline void SetPrintModulo(G4int value) {
        fPrintModulo = value;
    }
    /**
      * Add energy deposition and sum
      * \param edep energy in step in keV
      */
    inline void     AddEdep(G4double edep) {
        fEnergyPerEvent += edep;
    }
    /**
     * Get energy of the current event
     */
    inline G4double GetEdep() const {
        return fEnergyPerEvent;
    }
private:
    /**
     * Hits collection
     * \param hcID the hit ID
     * \param event the current event
     */
    DetectorHitsCollection *GetHitsCollection(G4int hcID, 
                                              const G4Event *event) const;
    /**
     * Prints primary event statistics
     * \param absoEdep absolute energy deposition
     * \param absoTrackLength absolute track length of primary track
     */
    void PrintEventStatistics(G4double absoEdep,
                              G4double absoTrackLength) const;

    /**
     * 
     * 
     */
    G4int  fSensorHCID; //FIXME can be removed !?
    G4int  fPrintModulo;
    // the name of the digitizer
    G4String digitizerName;
    // energy per event
    G4double fEnergyPerEvent;
    G4int count;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


