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
/// \file include/Messenger.hh
/// \brief Definition of a general messenger class

#include "globals.hh"
#include "G4AutoLock.hh"

class G4GenericMessenger;
/**
 * implementation of a Messenger class
 */
class PrimaryGeneratorMessenger
{
public:
   static PrimaryGeneratorMessenger *GetInstance();

    /**
     * Constructor of the Messenger class
     */
    PrimaryGeneratorMessenger();
    /**
     * Destructor of the Messenger class
     */
    virtual ~PrimaryGeneratorMessenger();
    /**
     * Sets the file name to
     * \param name file name string
     */
    void beamOn(G4int);
//     void setNumberOfFrames(G4int);
//     void SetNumberOfParticles(G4int); 	//for dacscan
//     void SetStepSize(G4double);		//for dacscan
    /**
     * Define the range for the dac scan in keV
     * \param start start dac value in keV
     * \param stop stop dac value in keV
     */
    void Dacscan(G4double, G4double);

    inline G4double GetStepSize() {
        return step;
    }
    inline G4int GetNumberOfParticles() {
        return nPar;
    }
    inline void setNumberOfFrames(G4int nf) {
        nFrames = nf;
    }

    inline G4int GetNumberOfFrames(){
        return nFrames; 
    }
    inline void SetStepSize(G4double sz) {
        step = sz;
    }

    inline void SetNumberOfParticles(G4int np) {
        nPar = np;
    }

private:
    /**
     * The messenger instance
     */
    static PrimaryGeneratorMessenger *instance;
    G4GenericMessenger *fMessenger;
    G4int       nFrames; /**< number for frames*/
    G4double    step; /**< step for dacscan */
    G4int       nPar; /**< number of particles for each step in the dacscan */


    /**
     * Initializes the messenger commands
     */
    void DefineCommands();
};
