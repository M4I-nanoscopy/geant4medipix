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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class

#include "DetectorMessenger.hh"
#include "DetectorConstructionBase.hh"
#include "MpxDetector.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "boost/filesystem.hpp"


DetectorMessenger::DetectorMessenger(DetectorConstructionBase *Det)
    : Detector(Det)
{
    pDir = new G4UIdirectory("/Medipix/");
    pDir->SetGuidance("UI command for Medipix chip and sensor");

    oldPath = "";
    
    G4String matList;
    const G4MaterialTable *matTbl = G4Material::GetMaterialTable();
    for (size_t i = 0; i < G4Material::GetNumberOfMaterials(); i++) {
        matList += (*matTbl)[i]->GetName();
        matList += " ";
    }
    matList += "G4_GALLIUM_ARSENIDE";
    matList += " ";

    //---------------------------------------------------------------------------- Sensor
    pSensorMatCmd = new G4UIcmdWithAString("/Sensor/material", this);
    pSensorMatCmd->SetGuidance("Select material of the Sensor");
    pSensorMatCmd->SetDefaultValue("G4_Si");
    pSensorMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    pSensorMatCmd->SetCandidates(matList);

    pSensorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/Sensor/thickness", this);
    pSensorThicknessCmd->SetGuidance("Set thickness of the Sensor");
    pSensorThicknessCmd->SetGuidance("Thickness in Z direction");
    pSensorThicknessCmd->SetDefaultValue(300.);
    pSensorThicknessCmd->SetDefaultUnit("um");
    pSensorThicknessCmd->SetParameterName("thickness", false);
    pSensorThicknessCmd->SetUnitCategory("Length");
    pSensorThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pNbPixelsCmd = new G4UIcmdWithAnInteger("/Sensor/pixels", this);
    pNbPixelsCmd->SetGuidance("Set number of pixels in each direction.");
    pNbPixelsCmd->SetParameterName("NbPixels", false);
    pNbPixelsCmd->SetDefaultValue(10);
    pNbPixelsCmd->SetRange("NbPixels>0");
    pNbPixelsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pPixelSizeCmd = new G4UIcmdWithADoubleAndUnit("/Sensor/pixelSize", this);
    pPixelSizeCmd->SetGuidance("Set pixel size with unit");
    pPixelSizeCmd->SetDefaultValue(110.);
    pPixelSizeCmd->SetDefaultUnit("um");
    pPixelSizeCmd->SetParameterName("PixelSize", false);
    pPixelSizeCmd->SetRange("PixelSize>0");
    pPixelSizeCmd->SetUnitCategory("Length");
    pPixelSizeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pRotationCmd = new G4UIcmdWithADouble("/Sensor/rotation", this);
    pRotationCmd->SetGuidance("Set the angle of rotation around the x-axis in degree");
    pRotationCmd->SetDefaultValue(0.);
    pRotationCmd->SetParameterName("Rotation", false);
    pRotationCmd->SetRange("0<=Rotation<=90");
    pRotationCmd->AvailableForStates(G4State_PreInit, G4State_Init);

    //---------------------------------------------------------------------------- Detector
    pSetDigitizerCmd = new G4UIcmdWithAString("/Detector/digitizer", this);
    pSetDigitizerCmd->SetGuidance("Select one Digitizer");
    pSetDigitizerCmd->SetDefaultValue("MpxDigitizer");
    pSetDigitizerCmd->SetParameterName("choice", false);
    pSetDigitizerCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    pSetDigitizerCmd->SetCandidates("MpxDigitizer WFDigitizer");

    pSetDetectorTypeCmd = new G4UIcmdWithAString("/Detector/type", this);
    pSetDetectorTypeCmd->SetGuidance("Select one Detector");
    pSetDetectorTypeCmd->SetDefaultValue("Medipix3RX");
    pSetDetectorTypeCmd->SetParameterName("choice", false);
    pSetDetectorTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    pSetDetectorTypeCmd->SetCandidates("Medipix3RX Timepix Timepix3 Dosepix");

    pSetDetectorGainCmd = new G4UIcmdWithAString("/Detector/gain", this);
    pSetDetectorGainCmd->SetGuidance("Select one capacity");
    pSetDetectorGainCmd->SetParameterName("choice", false);
    pSetDetectorGainCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    pSetDetectorGainCmd->SetCandidates("SLGM LGM HGM SHGM");

    pSetCsmModeCmd = new G4UIcmdWithABool("/Detector/csm", this);
    pSetCsmModeCmd->SetGuidance("Enable or disable the CSM Mode");
    pSetCsmModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    //---------------------------------------------------------------------------- Chip specific
    pSetThreshold1Cmd = new G4UIcmdWithADoubleAndUnit("/Medipix/threshold1", this);
    pSetThreshold1Cmd->SetGuidance("Set low threshold with unit");
    pSetThreshold1Cmd->SetParameterName("threshold1", false);
    pSetThreshold1Cmd->SetDefaultValue(4);
    pSetThreshold1Cmd->SetDefaultUnit("keV");
    pSetThreshold1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetThreshold2Cmd = new G4UIcmdWithADoubleAndUnit("/Medipix/threshold2", this);
    pSetThreshold2Cmd->SetGuidance("Set low threshold with unit");
    pSetThreshold2Cmd->SetParameterName("threshold2", false);
    pSetThreshold2Cmd->SetDefaultValue(4);
    pSetThreshold2Cmd->SetDefaultUnit("keV");
    pSetThreshold2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetThreshold3Cmd = new G4UIcmdWithADoubleAndUnit("/Medipix/threshold3", this);
    pSetThreshold3Cmd->SetGuidance("Set low threshold with unit");
    pSetThreshold3Cmd->SetParameterName("threshold3", false);
    pSetThreshold3Cmd->SetDefaultValue(4);
    pSetThreshold3Cmd->SetDefaultUnit("keV");
    pSetThreshold3Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetThreshold4Cmd = new G4UIcmdWithADoubleAndUnit("/Medipix/threshold4", this);
    pSetThreshold4Cmd->SetGuidance("Set low threshold with unit");
    pSetThreshold4Cmd->SetParameterName("threshold4", false);
    pSetThreshold4Cmd->SetDefaultValue(4);
    pSetThreshold4Cmd->SetDefaultUnit("keV");
    pSetThreshold4Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetTpxThresholdCmd = new G4UIcmdWithADoubleAndUnit("/Timepix/threshold", this);
    pSetTpxThresholdCmd->SetGuidance("Set low threshold with unit");
    pSetTpxThresholdCmd->SetParameterName("threshold", false);
    pSetTpxThresholdCmd->SetDefaultValue(4);
    pSetTpxThresholdCmd->SetDefaultUnit("keV");
    pSetTpxThresholdCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetTpxModeCmd = new G4UIcmdWithAString("/Timepix/mode", this);
    pSetTpxModeCmd->SetGuidance("Set operation mode of Timepix");
    pSetTpxModeCmd->SetParameterName("tpxMode", false);
    pSetTpxModeCmd->SetDefaultValue("PC");
    pSetTpxModeCmd->SetCandidates("PC ToT ToA");
    pSetTpxModeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pUpdateCmd = new G4UIcmdWithoutParameter("/Sensor/Update", this);
    pUpdateCmd->SetGuidance("Update geometry.");
    pUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
    pUpdateCmd->SetGuidance("if you changed geometrical value(s)");
    pUpdateCmd->AvailableForStates(G4State_PreInit, G4State_Idle);



    //---------------------------------------------------------------------------- Filter
    pSetFilterCmd = new G4UIcmdWithABool("/Filter/use", this);
    pSetFilterCmd->SetGuidance("Place filter in front of sensor");
    pSetFilterCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pFilterMatCmd = new G4UIcmdWithAString("/Filter/material", this);
    pFilterMatCmd->SetGuidance("Select material of the filter");
    pFilterMatCmd->SetDefaultValue("G4_Al");
    pFilterMatCmd->SetParameterName("choice", false);
    pFilterMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    pFilterMatCmd->SetCandidates(matList);

    pFilterThicknessCmd = new G4UIcmdWithADoubleAndUnit("/Filter/thickness", this);
    pFilterThicknessCmd->SetGuidance("Set thickness of the filter");
    pFilterThicknessCmd->SetGuidance("Thickness in Z direction");
    pFilterThicknessCmd->SetDefaultValue(300.);
    pFilterThicknessCmd->SetDefaultUnit("um");
    pFilterThicknessCmd->SetParameterName("thickness", false);
    pFilterThicknessCmd->SetUnitCategory("Length");
    pFilterThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pFilterZCmd = new G4UIcmdWithADoubleAndUnit("/Filter/z", this);
    pFilterZCmd->SetGuidance("Set thickness z position of the filter edge relative to the sensor");
    pFilterZCmd->SetDefaultValue(300.);
    pFilterZCmd->SetDefaultUnit("um");
    pFilterZCmd->SetParameterName("z", false);
    pFilterZCmd->SetUnitCategory("Length");
    pFilterZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pFilterRotationCmd = new G4UIcmdWithADouble("/Filter/rotation", this);
    pFilterRotationCmd->SetGuidance("Set Rotation of the filter around the x-axis in degree.");
    pFilterRotationCmd->SetDefaultValue(0.);
    pFilterRotationCmd->SetParameterName("Rotation", false);
    pFilterRotationCmd->SetRange("0<=Rotation<=90");
    pFilterRotationCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    //---------------------------------------------------------------------------- Bump bonds
    pSetBumpsCmd = new G4UIcmdWithABool("/Bumps/use", this);
    pSetBumpsCmd->SetGuidance("Turn on or off bump bonds");
    pSetBumpsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pBumpRadiiCmd = new G4UIcmdWithADoubleAndUnit("/Bumps/radii", this);
    pBumpRadiiCmd->SetGuidance("Set the radii of the bumps");
    pBumpRadiiCmd->SetDefaultValue(30.);
    pBumpRadiiCmd->SetDefaultUnit("um");
    pBumpRadiiCmd->SetParameterName("r", false);
    pBumpRadiiCmd->SetUnitCategory("Length");
    pBumpRadiiCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pBumpHeightCmd = new G4UIcmdWithADoubleAndUnit("/Bumps/height", this);
    pBumpHeightCmd->SetGuidance("Set the height of the bumps");
    pBumpHeightCmd->SetDefaultValue(50.);
    pBumpHeightCmd->SetDefaultUnit("um");
    pBumpHeightCmd->SetParameterName("h", false);
    pBumpHeightCmd->SetUnitCategory("Length");
    pBumpHeightCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetBumpMaterialCmd = new G4UIcmdWithAString("/Bumps/material", this);
    pSetBumpMaterialCmd->SetGuidance("Select material of the bump bonds");
    pSetBumpMaterialCmd->SetDefaultValue("G4_Cu");
    pSetBumpMaterialCmd->SetParameterName("choice", false);
    pSetBumpMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    pSetBumpMaterialCmd->SetCandidates(matList);


    //---------------------------------------------------------------------------- Additional Solids
    pSetChipCmd = new G4UIcmdWithABool("/Solids/chip", this);
    pSetChipCmd->SetGuidance("Place electronics chip behind sensor and bumps");
    pSetChipCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetCollimatorCmd = new G4UIcmdWithABool("/Solids/Collimator/use", this);
    pSetCollimatorCmd->SetGuidance("Place collimator in front of each pixel");
    pSetCollimatorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetPcbCmd = new G4UIcmdWithABool("/Solids/pcb", this);
    pSetPcbCmd->SetGuidance("Place pcb");
    pSetPcbCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetWorldMatCmd = new G4UIcmdWithAString("/Solids/world", this);
    pSetWorldMatCmd->SetGuidance("Select material of the world");
    pSetWorldMatCmd->SetDefaultValue("G4_Galactic");
    pSetWorldMatCmd->SetParameterName("choice", false);
    pSetWorldMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    pSetWorldMatCmd->SetCandidates(matList);

    pSetCollimatordMatCmd = new G4UIcmdWithAString("/Solids/Collimator/material", this);
    pSetCollimatordMatCmd->SetGuidance("Select material of the collimator");
    pSetCollimatordMatCmd->SetDefaultValue("G4_Fe");
    pSetCollimatordMatCmd->SetParameterName("choice", false);
    pSetCollimatordMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    pSetCollimatordMatCmd->SetCandidates(matList);


    pSetCollimatorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/Solids/Collimator/thickness", this);
    pSetCollimatorThicknessCmd->SetGuidance("Set the thickness of the collimator");
    pSetCollimatorThicknessCmd->SetDefaultValue(500.);
    pSetCollimatorThicknessCmd->SetDefaultUnit("um");
    pSetCollimatorThicknessCmd->SetParameterName("h", false);
    pSetCollimatorThicknessCmd->SetUnitCategory("Length");
    pSetCollimatorThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetCollimatorRadiiCmd = new G4UIcmdWithADoubleAndUnit("/Solids/Collimator/radii", this);
    pSetCollimatorRadiiCmd->SetGuidance("Set the radii of the collimator holes");
    pSetCollimatorRadiiCmd->SetDefaultValue(55.);
    pSetCollimatorRadiiCmd->SetDefaultUnit("um");
    pSetCollimatorRadiiCmd->SetParameterName("h", false);
    pSetCollimatorRadiiCmd->SetUnitCategory("Length");
    pSetCollimatorRadiiCmd->AvailableForStates(G4State_PreInit, G4State_Idle);



    pSetBlockCmd = new G4UIcmdWithABool("/Solids/block", this);
    pSetBlockCmd->SetGuidance("Place a metal block");
    pSetBlockCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetBlockZPosCmd = new G4UIcmdWithADoubleAndUnit("/Solids/block/z", this);
    pSetBlockZPosCmd->SetGuidance("Set the distance between the chip and the radioactive block");
    pSetBlockZPosCmd->SetDefaultValue(55.);
    pSetBlockZPosCmd->SetDefaultUnit("um");
    pSetBlockZPosCmd->SetParameterName("h", false);
    pSetBlockZPosCmd->SetUnitCategory("Length");
    pSetBlockZPosCmd->AvailableForStates(G4State_PreInit, G4State_Idle);


    pSetBlockThicknessCmd = new G4UIcmdWithADoubleAndUnit("/Solids/block/thickness", this);
    pSetBlockThicknessCmd->SetGuidance("Set the distance between the chip and the radioactive block");
    pSetBlockThicknessCmd->SetDefaultValue(55.);
    pSetBlockThicknessCmd->SetDefaultUnit("um");
    pSetBlockThicknessCmd->SetParameterName("h", false);
    pSetBlockThicknessCmd->SetUnitCategory("Length");
    pSetBlockThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    
    pSetBlockWidthCmd = new G4UIcmdWithADoubleAndUnit("/Solids/block/width", this);
    pSetBlockWidthCmd->SetGuidance("Set the width of the radioactive block");
    pSetBlockWidthCmd->SetDefaultValue(55.);
    pSetBlockWidthCmd->SetDefaultUnit("um");
    pSetBlockWidthCmd->SetParameterName("h", false);
    pSetBlockWidthCmd->SetUnitCategory("Length");
    pSetBlockWidthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    //----------------------------------------------------------------------------Output
    pSetOutputFilenameCmd = new G4UIcmdWithAString("/Output/frames", this);
    pSetOutputFilenameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetSparseOutputFilenameCmd = new G4UIcmdWithAString("/Output/sparse", this);
    pSetSparseOutputFilenameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetOutputDirectoryCmd = new G4UIcmdWithAString("/Output/path", this);
    pSetOutputDirectoryCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetConfigFilenameCmd = new G4UIcmdWithAString("/Output/config", this);
    pSetConfigFilenameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    pSetStoreTrajCmd = new G4UIcmdWithABool("/Output/store", this);
    pSetStoreTrajCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{

    //Sensor
    delete pSensorMatCmd;
    delete pSensorThicknessCmd;
    delete pNbPixelsCmd;
    delete pPixelSizeCmd;
    delete pRotationCmd;

    //Detector
    delete pSetDigitizerCmd;
    delete pSetDetectorTypeCmd;
    delete pSetDetectorGainCmd;
    //delete pSetThDispCmd;
    delete pSetCsmModeCmd;
    //delete pSetNoiseCmd;

    //Other
    delete pUpdateCmd;
    delete pDir;

    //Chip specific
    delete pSetThreshold1Cmd;
    delete pSetThreshold2Cmd;
    delete pSetThreshold3Cmd;
    delete pSetThreshold4Cmd;
    delete pSetTpxThresholdCmd;

    //Filters
    delete pSetFilterCmd;
    delete pFilterMatCmd;
    delete pFilterThicknessCmd;
    delete pFilterZCmd;
    delete pFilterRotationCmd;

    //Solids
    delete pSetCollimatorCmd;
    delete pSetCollimatordMatCmd;
    delete pSetCollimatorThicknessCmd;
    delete pSetCollimatorRadiiCmd;
    delete pSetChipCmd;
    delete pSetPcbCmd;
    delete pSetWorldMatCmd;
    delete pSetBlockCmd;
    delete pSetBlockZPosCmd;
    delete pSetBlockThicknessCmd;
    delete pSetBlockWidthCmd;


    //Bumps
    delete pSetBumpsCmd;
    delete pBumpRadiiCmd;
    delete pBumpHeightCmd;

    //Output
    delete pSetOutputFilenameCmd;
    delete pSetSparseOutputFilenameCmd;
    delete pSetOutputDirectoryCmd;
    delete pSetConfigFilenameCmd;
    delete pSetStoreTrajCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand *cmd, G4String newValue)
{

    // Sensor
    if (cmd == pSensorMatCmd) {
        Detector->SetSensorMaterial(newValue);
    }
    if (cmd == pSensorThicknessCmd) {
        Detector->SetSensorThickness(pSensorThicknessCmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pNbPixelsCmd) {
        Detector->SetNbPixels(pNbPixelsCmd->GetNewIntValue(newValue));

        MpxDetector *mpxDet = MpxDetector::GetInstance();
        mpxDet->ReInitMatrix();
    }
    if (cmd == pPixelSizeCmd) {
        Detector->SetPixelSize(pPixelSizeCmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pRotationCmd) {
        Detector->SetRotation(pRotationCmd->GetNewDoubleValue(newValue));
    }

    // Detector
    if (cmd == pSetDigitizerCmd) {
        Detector->SetDigitizerName(newValue);
    }
    if (cmd == pSetDetectorTypeCmd) {
        Detector->SetDetectorType(newValue);
    }
    if (cmd == pSetDetectorGainCmd) {
        Detector->SetDetectorGain(newValue);
    }
    if (cmd == pSetCsmModeCmd) {
        Detector->SetCsmMode(pSetCsmModeCmd->GetNewBoolValue(newValue));
    }

    // Chip specific
    if (cmd == pSetThreshold1Cmd) {
        MpxDetector *mpxDet = MpxDetector::GetInstance();
        mpxDet->SetMpxThreshold1(pSetThreshold1Cmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pSetThreshold2Cmd) {
        MpxDetector *mpxDet = MpxDetector::GetInstance();
        mpxDet->SetMpxThreshold2(pSetThreshold2Cmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pSetThreshold3Cmd) {
        MpxDetector *mpxDet = MpxDetector::GetInstance();
        mpxDet->SetMpxThreshold3(pSetThreshold3Cmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pSetThreshold4Cmd) {
        MpxDetector *mpxDet = MpxDetector::GetInstance();
        mpxDet->SetMpxThreshold4(pSetThreshold4Cmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pSetTpxThresholdCmd) {
        MpxDetector *mpxDet = MpxDetector::GetInstance();
        mpxDet->SetTpxThreshold(pSetTpxThresholdCmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pSetTpxModeCmd) {
        Detector->SetTpxMode(newValue);
    }

    //Other
    if (cmd == pUpdateCmd) {
        Detector->UpdateGeometry();
    }

    //Filters
    if (cmd == pSetFilterCmd) {
        Detector->SetFilterBool(pSetFilterCmd->GetNewBoolValue(newValue));
    }
    if (cmd == pFilterMatCmd) {
        Detector->SetFilterMaterial(newValue);
    }
    if (cmd == pFilterThicknessCmd) {
        Detector->SetFilterThickness(pFilterThicknessCmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pFilterZCmd) {
        Detector->SetFilterZ(pFilterZCmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pFilterRotationCmd){
        Detector->SetFilterRotation(pFilterRotationCmd->GetNewDoubleValue(newValue));
    }


    //Bumps
    if (cmd == pSetBumpsCmd) {
        Detector->SetBumpsBool(pSetBumpsCmd->GetNewBoolValue(newValue));
    }
    if (cmd == pBumpRadiiCmd) {
        Detector->SetBumpRadii(pBumpRadiiCmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pBumpHeightCmd) {
        Detector->SetBumpHeight(pBumpHeightCmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pSetBumpMaterialCmd) {
        Detector->SetBumpMaterial(newValue);
    }

    //Solids
    if (cmd == pSetCollimatorCmd) {
        Detector->SetCollimatorBool(pSetCollimatorCmd->GetNewBoolValue(newValue));
    }
    if (cmd == pSetCollimatordMatCmd) {
        Detector->SetCollimatorMaterial(newValue);
    }
    if (cmd == pSetCollimatorThicknessCmd) {
        Detector->SetCollimatorThickness(pSetCollimatorThicknessCmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pSetCollimatorRadiiCmd) {
        Detector->SetCollimatorRadii(pSetCollimatorRadiiCmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pSetChipCmd) {
        Detector->SetChipBool(pSetChipCmd->GetNewBoolValue(newValue));
    }
    if (cmd == pSetPcbCmd) {
        Detector->SetPcbBool(pSetPcbCmd->GetNewBoolValue(newValue));
    }
    if (cmd == pSetWorldMatCmd) {
        Detector->SetWorldMaterial(newValue);
    }
    if (cmd == pSetBlockCmd) {
        Detector->SetBlockBool(pSetBlockCmd->GetNewBoolValue(newValue));
    }
    if (cmd == pSetBlockZPosCmd) {
        Detector->SetBlockZ(pSetBlockZPosCmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pSetBlockThicknessCmd) {
        Detector->SetBlockThickness(pSetBlockThicknessCmd->GetNewDoubleValue(newValue));
    }
    if (cmd == pSetBlockWidthCmd) {
        Detector->SetBlockWidth(pSetBlockWidthCmd->GetNewDoubleValue(newValue));
    }

    //Output
    if (cmd == pSetOutputFilenameCmd) {
        Detector->SetOutputFilename(newValue);
    }
    if (cmd == pSetSparseOutputFilenameCmd) {
        Detector->SetSparseOutputFilename(newValue);
    }
    if (cmd == pSetOutputDirectoryCmd) {
        // change directory and try to create it if not there
        G4cout << "Trying to change to: " << newValue << G4endl;
        try {
	  if (oldPath == ""){
	      oldPath = boost::filesystem::current_path().string(); //save current wd
	  }else{
	    boost::filesystem::current_path( oldPath );
	  }
        boost::filesystem::current_path(newValue);
        } catch (const boost::filesystem::filesystem_error &e) {
            if (e.code() == boost::system::errc::no_such_file_or_directory) {
                G4cout << "Directory does not exist. Trying to create." << G4endl;
                try {
                    boost::filesystem::create_directories(newValue);
                    boost::filesystem::current_path(newValue);
                    G4cout << "Working directory is now: " << boost::filesystem::current_path() << G4endl;
                } catch (const boost::filesystem::filesystem_error &e2) {
                    G4cout << e2.code().message() << G4endl;
                    G4cout << "Working directory remains: " << boost::filesystem::current_path() << G4endl;
                }
            }
        }
    }
    if (cmd == pSetConfigFilenameCmd) {
        Detector->SetConfigFilename(newValue);
    }
    if (cmd == pSetStoreTrajCmd) {
        Detector->SetStoreTraj(pSetStoreTrajCmd->GetNewBoolValue(newValue));
    }
}
