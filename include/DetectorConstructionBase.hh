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
/// \file include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstructionBase class

#ifndef DetectorConstructionBase_h
#define DetectorConstructionBase_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <G4RotationMatrix.hh>

class G4Box;
class G4Material;
class G4VPhysicalVolume;
class G4GenericMessenger;
class DetectorMessenger;

class DetectorConstructionBase : public G4VUserDetectorConstruction
{
public:
    DetectorConstructionBase();
    virtual ~DetectorConstructionBase();

public:
    virtual G4VPhysicalVolume *Construct();
    virtual void ConstructSDandField();

    void UpdateGeometry();

    // set methods
    //
    void SetWorldMaterial(const G4String &);
    void SetSensorMaterial(const G4String &);
    void SetFilterMaterial(const G4String &);
    void SetFilterRotation(G4double xangle);
    void SetSensorThickness(G4double);
    void SetNbPixels(G4int);
    void SetPixelSize(G4double);
    void SetTpxMode(G4String);
    inline void SetDigitizerName(G4String digitizer) {
        digitizerName = digitizer;
    }
    void SetDetectorType(G4String);
    void SetDetectorGain(G4String);
    void SetRotation(G4double);
    inline void SetCsmMode(G4bool mode) {
        csmMode = mode;
    }

    inline void SetNoise(G4double noise) {
        eNoise = noise;
    }

    //Filter
    /**
     * 
     */
    inline void SetFilterBool(G4double fi) {
        filter = fi;
    }
    // filter thickness, fluorescence plate
    void SetFilterThickness(G4double);
    // position of the filter/fluorescence plate
    void SetFilterZ(G4double);

    //Bumps
    // enable or disble (default)
    inline void SetBumpsBool(G4double bm) {
        bumps = bm;
    }
    // bump radius (cylinder)
    inline void SetBumpRadii(G4double r) {
        bumpRadii = r;
    }
    // bump height
    inline void SetBumpHeight(G4double h) {
        bumpHeight = h;
    }
    // material of the bump
    void SetBumpMaterial(const G4String &);

    // Collimators
    inline void SetCollimatorBool(G4double cl) {
        collimator = cl;
    }
    void SetCollimatorMaterial(const G4String &);

    inline void SetCollimatorThickness(G4double h) {
        collimatorThickness = h;
    }
    inline void SetCollimatorRadii(G4double h) {
        colRadii = h;
    }
    // extra solids
    // add the chip; fixed material, fixed dimensions
    inline void SetChipBool(G4bool ch) {
        chip = ch;
    }
    // add the PCB
    inline void SetPcbBool(G4bool ch) {
        pcb = ch;
    }

    //Block FIXME: what is that
    inline void SetBlockBool(G4bool ch) {
        block = ch;
    }
    inline void SetBlockZ(G4double h) {
        actPosZ = h;
    }
    inline void SetBlockThickness(G4double h) {
        actZ = h;
    }
    inline void SetBlockWidth(G4double h) {
        actX = h;
	actY = h;
    }

    //Output
    inline void SetConfigFilename(G4String f) {
        configFilename = f;
    }
    inline void SetOutputFilename(G4String fname) {
        optFname = fname;
    }
    inline void SetSparseOutputFilename(G4String fname) {
        sparseOptFname = fname;
    }

    //
    virtual G4VPhysicalVolume* DefineVolumes() = 0;
    void CalculateGeometry();


    G4int       GetNbPixels() {
        return nPixel;
    }
    G4double    GetPixelSize() {
        return pixelSize;
    }
    G4double    GetSensorThickness() {
        return sensorThickness;
    }
    G4Material *GetSensorMaterial() {
        return sensorMaterial;
    }
    G4Material *GetFilterMaterial() {
        return filterMaterial;
    }
    G4int GetTpxMode()	{
        return 	tpxMode;
    }
    G4String GetDigitizerName() {
        return digitizerName;
    };
    G4int GetDetectorType() {
        return detectorType;
    };
    G4double GetDetectorGain() {
        return feedbackCapacitance;
    };

    G4bool GetCsmMode() {
        return csmMode;
    }
    G4double GetNoise() {
        return eNoise;
    }
    G4double GetFilterBool() {
        return filter;
    }
    G4double GetFilterThickness() {
        return filterThickness;
    }
    G4double GetFilterZ() {
        return filterZ;
    }

    //----------------------------------Bump bonds
    G4double GetBumpsBool() {
        return bumps;
    }
    G4double GetBumpRadii() {
        return bumpRadii;
    }
    G4double GetBumpHeight() {
        return bumpHeight;
    }
    G4Material *GetBumpMaterial() {
        return bumpMaterial;
    }

    //-------------------------------------Solids
    G4double GetCollimatorBool() {
        return collimator;
    }
    G4Material *GetCollimatorMaterial() {
        return collimatorMaterial;
    }
    G4double GetCollimatorThickness() {
        return collimatorThickness;
    }
    G4double GetCollimatorRadii() {
        return colRadii;
    }
    G4double GetChipBool() {
        return chip;
    }
    G4double GetPcbBool() {
        return pcb;
    }
    G4Material *GetWorldMaterial() {
        return worldMaterial;
    }
    G4bool GetBlockBool() {
        return block;
    }
    G4double GetBlockZ() {
        return actPosZ;
    }
    G4double GetBlockThickness() {
        return actZ;
    }
    //Output
    G4String GetOutputFilename() {
        return optFname;
    }
    G4String GetSparseOutputFilename() {
        return sparseOptFname;
    }
    G4String GetConfigFilename() {
        return configFilename;
    }


    static G4ThreadLocal G4GenericMessenger *fMessenger;  // messenger

    G4double    sensorSizeXY;
    G4bool      fCheckOverlaps; // option to activate checking of volumes overlaps
    G4Material *sensorMaterial;

    G4RotationMatrix fRotation;
    G4double    sensorThickness;
    G4int       nPixel;
    G4double    pixelSize;

    //Filters
    G4bool 	filter;
    G4Material  *filterMaterial;
    G4double	filterThickness;
    G4double 	filterZ;
    G4RotationMatrix filterRotation;

    //Bump bonds
    G4bool 	bumps;
    G4Material *bumpMaterial;
    G4double bumpRadii;
    G4double bumpHeight;

    //Solids
    G4bool	collimator;
    G4Material *collimatorMaterial;
    G4double collimatorThickness;
    G4double colRadii;
    G4bool chip;
    G4bool pcb;
    G4Material *worldMaterial;
    G4double block;
    G4double actX;
    G4double actY;
    G4double actZ;
    G4double actPosZ;

private:
    //
    G4String    digitizerName;
    G4double    feedbackCapacitance;
    G4double    eNoise;
    G4int       detectorType;
    G4int       tpxMode;
    G4bool      csmMode;

    //Output
    G4String    optFname;
    G4String	sparseOptFname;
    G4String    configFilename;

    DetectorMessenger *fDetectorMessenger;
    void DefineMaterials();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

