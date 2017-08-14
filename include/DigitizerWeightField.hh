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
/// \file include/WeightFieldDigitizer.hh
/// \brief Definition of the WeightFieldDigitizer class

#ifndef DigitizerWeightField_h
#define DigitizerWeightField_h 1

#include "DetectorHit.hh"
#include "MpxDetector.hh"

#include "G4VDigitizerModule.hh"
#include "ExportMgr.hh"
#include <vector>
#include <boost/concept_check.hpp>
#include "PreampMedipix.hh"

class DetectorConstructionBase;

class DigitizerWeightField : public G4VDigitizerModule
{
public:
    DigitizerWeightField(G4String aName);

    ~DigitizerWeightField();

    virtual void Digitize();
 

private:
    MpxDigitCollection *digitCollection;

    /**
     * sensor thickness
     */
    G4double detectorThickness;
    /**
     * reverse bias voltage of the sensor
     */ 
    G4double biasVoltage;
    /**
     * depletion voltage for overdepleted sensor
     */
    G4double depletionVoltage;
    /**
     * sensor depletion length
     */
    G4double depletedDepth;
    /**
     * experiment temperature (influences mobilities)
     */
    G4double Temperature;
    /**
     * number of pixels in x direction
     */
    G4int    nPixX;
    /**
     * number of pixels in y direction
     */
    G4int    nPixY;
    /**
     * size of the pixels in um
     */
    G4double pixelSize;
    /**
     * number of electron hole pairs from energy deposition
     */
    G4double nElectronHolePairs;
    
    G4double Default_Relative_Permittivity;
    /**
     * chose between electron 0, and hole collection 1
     */
    G4bool   typeToCollect;
    /**
     * collect electrons and holes if true
     */
    G4bool   trackBothTypes;
    /**
     * amount of charge to track in one cloud
     */
    G4int    nChargeToTrackTogether;
    /**
     * place subcharges in sphere around the interaction?
     */
    G4bool   initialDisplacement;
    /**
     * charge cloud size
     */
    G4double chargeCloudSigma;
    /**
     * amount of electronics noise in elementary charges
     */
    G4double elecSigma;
    /**
     * fano factor
     */
    G4double fanoFactor;
    /**
     * use the Diffusion/Repulsion model?
     */
    G4bool   useDiffusionRepulsion;
    /**
     * apply trapping?
     */
    G4bool   doTrapping;
    /**
     * trapping time
     */
    G4double trappingTime;
    /**
     * minimal energy in a pixel for the algorithm to cut
     */
    G4double minEnergy;
    /**
     * calculated mobility: depends on temperature...
     */
    G4double mobility ;
    /**
     * electric field in 3D
     */
    G4double electricFieldX;
    G4double electricFieldY;
    G4double electricFieldZ;
    
    // Default mobilities
    G4double Default_Electron_Mobility;
    G4double Default_Hole_Mobility;
    G4double Default_Electron_D;
    G4double Default_Hole_D;
    G4double Default_LET;
    


    //mobility dependence on electric field
    G4double Electron_AlphaField;
    G4double Electron_ThetaField;
    G4double Electron_TempNominal;
    G4double Electron_Beta;
    G4double Electron_Saturation_Velocity;
//     Electron_Saturation_Velocity = Electron_AlphaField*TMath::Power(1.0+Electron_ThetaField*TMath::Exp(Temperature/Electron_TempNominal),-1.0);
    G4double Hole_AlphaField;
    G4double Hole_ThetaField;
    G4double Hole_TempNominal;
    G4double Hole_Beta;
    G4double Hole_Saturation_Velocity;
    
    G4double pulsePrecision;
    G4double maxPulseTime;
    G4int    nPulseArrayElements;
    G4double ampResponseTime;
    G4int    nAmpResponseElements;
    G4double cutDepthTracking;
    G4bool   writePeakToFile;
    G4bool   useCSM;

    MpxDetector* detector;
    DetectorConstructionBase* myDet;
    
    G4int nGridPointsXY;
    G4int nGridPointsZ;
    G4double* potentialFine;
    G4double gridSizeFineXY;
    G4double gridSizeFineZ;
    /**
     * Creates weighting potential from file and detector geometry
     */
    void        CreatePotentialTable();
    G4double    GetInterpolatedPotential(G4ThreeVector);
    G4double    GetPotential(G4int x, G4int y, G4int z);
    /**
     * Compute the electric field in 3 dimensions
     * \param x coordinate
     * \param y coordinate
     * \param z coordinate
     */
    void ComputeElectricField(G4double /*x*/, G4double /*y*/, G4double /*z*/);
    /**
     * Compute the electric field in 1D, usually z
     * \param z the coordinate
     */
    void ComputeElectricField1D(G4double);
    /**
     * compute norm of the electic field
     * \param z coordinate
     */
    G4double GetElectricFieldNorm(G4double /*x*/, G4double /*y*/, G4double /*z*/);
    
    std::vector<G4double>  ComputeDrift(G4ThreeVector position, G4int column, G4int line, G4double time, G4double charge, std::map<std::pair<G4int,G4int>, G4double* >* pInducedPixelContent, G4bool particle, G4bool primType);
    
    std::vector<G4double>  RKF45Integration(G4double x, G4double y, G4double z,G4double dt, G4bool particle);
    
    G4double    MobilityElectronHole(G4double, G4double, G4double, G4bool particle);
    
    G4ThreeVector GetPixelCoordinates(G4ThreeVector pos);
    
    std::pair<G4int, G4int> GetPixelFromPosition(G4ThreeVector pos);
    /**
     * compute fano noise for given charge
     * \param charge amount of charge
     */
    G4double FanoNoise(G4double charge);
    /**
     * compute sigma of initial charge cloud
     * \param energy amount of energy of interaction step
     */
    void ComputeChargeCloudSigma(G4double energy);
    /**
     * compute diffusion repulsion sigma at tracking step
     * \param time tracking time
     * \param sigma last step sigma
     * \param charge amount of charge tracked
     * \param particle particletype
     */
    G4double ComputeDiffusionRepulsionSigma(G4double time, G4double sigma, G4double charge, G4bool particle);

    G4bool DEBUG;
    /**
     * @brief only charge collection, no weighting field.
     */
    G4bool chargeCollectionMode;

    PreampMedipix* preamp;
    
    G4String material;
    
    protected:
    //void fmod(double arg1, gridSizeXY arg2);
};
#endif
