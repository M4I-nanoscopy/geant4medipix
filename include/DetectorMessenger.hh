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
/// \file include/DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include <boost/concept_check.hpp>


class DetectorConstructionBase;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;

class DetectorMessenger: public G4UImessenger
{
public:
    DetectorMessenger(DetectorConstructionBase *);
    virtual ~DetectorMessenger();
    virtual void SetNewValue(G4UIcommand *, G4String);
    

private:
    DetectorConstructionBase        *Detector;
    G4UIdirectory               *pDir;
    G4UIcmdWithAString          *pSensorMatCmd;
    G4UIcmdWithAString          *pSetDigitizerCmd;
    G4UIcmdWithAString          *pSetDetectorTypeCmd;
    G4UIcmdWithAnInteger        *pNbPixelsCmd;
    G4UIcmdWithABool            *verboseCmd;
    G4UIcmdWithADoubleAndUnit   *pSensorThicknessCmd;
    G4UIcmdWithoutParameter     *pUpdateCmd;
    G4UIcmdWithADoubleAndUnit   *pPixelSizeCmd;
    G4UIcmdWithADoubleAndUnit   *pSetThreshold1Cmd;
    G4UIcmdWithADoubleAndUnit   *pSetThreshold2Cmd;
    G4UIcmdWithADoubleAndUnit   *pSetThreshold3Cmd;
    G4UIcmdWithADoubleAndUnit   *pSetThreshold4Cmd;
    G4UIcmdWithADoubleAndUnit   *pSetTpxThresholdCmd;
    G4UIcmdWithAString          *pSetTpxModeCmd;
    G4UIcmdWithADouble          *pRotationCmd;
    
    G4UIcmdWithAString *pSetDetectorGainCmd;
    // detector modes, noise
    G4UIcmdWithADouble  *pSetThDispCmd;
    G4UIcmdWithABool    *pSetCsmModeCmd;
    G4UIcmdWithADouble  *pSetNoiseCmd;
    // extra solids, bumps, filters
    G4UIcmdWithABool          *pSetFilterCmd;
    G4UIcmdWithAString        *pFilterMatCmd;
    G4UIcmdWithADoubleAndUnit *pFilterThicknessCmd;
    G4UIcmdWithADoubleAndUnit *pFilterZCmd;
    G4UIcmdWithADouble        *pFilterRotationCmd;
    
    G4UIcmdWithABool            *pSetBumpsCmd;
    G4UIcmdWithADoubleAndUnit   *pBumpRadiiCmd;
    G4UIcmdWithADoubleAndUnit   *pBumpHeightCmd;
    G4UIcmdWithAString          *pSetBumpMaterialCmd;
    
    G4UIcmdWithABool            *pSetCollimatorCmd;
    G4UIcmdWithAString          *pSetCollimatordMatCmd;
    G4UIcmdWithADoubleAndUnit   *pSetCollimatorThicknessCmd;
    G4UIcmdWithADoubleAndUnit   *pSetCollimatorRadiiCmd;
    G4UIcmdWithABool            *pSetChipCmd;
    G4UIcmdWithABool            *pSetPcbCmd;
    G4UIcmdWithAString          *pSetWorldMatCmd;
    G4UIcmdWithABool            *pSetBlockCmd;
    G4UIcmdWithADoubleAndUnit   *pSetBlockZPosCmd;
    G4UIcmdWithADoubleAndUnit   *pSetBlockThicknessCmd;
    G4UIcmdWithADoubleAndUnit   *pSetBlockWidthCmd;
    
    
    //Output
    G4UIcmdWithAString  *pSetOutputFilenameCmd;
    G4UIcmdWithAString  *pSetSparseOutputFilenameCmd;
    G4UIcmdWithAString  *pSetOutputDirectoryCmd;
    G4UIcmdWithAString  *pSetConfigFilenameCmd;
    G4UIcmdWithAString  *pSetHdf5FilenameCmd;
    G4UIcmdWithABool *pSetStoreTrajCmd;

    G4String 		oldPath;
    
};

#endif
