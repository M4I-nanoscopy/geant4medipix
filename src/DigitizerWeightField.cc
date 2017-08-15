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
// $Id: WeightFieldDigitizer.cc $
//
/// \file src/WeightFieldDigitizer.cc
/// \brief Implementation of the WeightFieldDigitizer class


#include "DigitizerWeightField.hh"
#include "PreampMedipix.hh"
#include "DetectorHit.hh"
#include "DetectorConstructionBase.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4DigiManager.hh"
#include "G4SystemOfUnits.hh"

#include "CLHEP/Random/RandGauss.h"


#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"


using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DigitizerWeightField::DigitizerWeightField(G4String aName) : G4VDigitizerModule(aName)
{
    // Make the collectionName available to the DigitManager
    G4String colName = "DigitsCollection";
    collectionName.push_back(colName);

    //DetectorConstruction* detectorConstruction;
    G4RunManager *fRM = G4RunManager::GetRunManager();
    myDet = (DetectorConstructionBase *)(fRM->GetUserDetectorConstruction());

    //detector
    detector = MpxDetector::GetInstance();
    nPixX = myDet->GetNbPixels();
    nPixY = nPixX;
    pixelSize = myDet->GetPixelSize();
    detectorThickness = myDet->GetSensorThickness();
    elecSigma =  myDet->GetNoise();

    // weighting potential
    CreatePotentialTable();

    //load .ini file with configuration data
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("DetectorConfig.ini", pt);


    // get material and load material properties from ini file
    material =  myDet->GetSensorMaterial()->GetName();
    G4String sensor;
    if (material == "G4_Si")
        sensor = "sensor_silicon";
    else if (material == "G4_CADMIUM_TELLURIDE")
    {
        sensor = "sensor_cdte";
    }
    //load sensor properties
    biasVoltage = pt.get<G4double>(sensor + ".biasVoltage");
    depletionVoltage = pt.get<G4double>(sensor + ".depletionVoltage");

    Temperature = pt.get<G4double>(sensor + ".Temperature");
    depletedDepth = detectorThickness; //TODO make it dynamic


    if (material == "G4_CADMIUM_TELLURIDE")
    {
      electricFieldZ = biasVoltage/detectorThickness;
    }
    //Silicon
    nElectronHolePairs = pt.get<G4double>(sensor + ".nElectronHolePairs");
    fanoFactor = pt.get<G4double>(sensor + ".fanoFactor");

    Default_Relative_Permittivity = pt.get<G4double>(sensor + ".Default_Relative_Permittivity");
    typeToCollect = pt.get<G4bool>("computation.typeToCollect");                          //0: electron || 1: hole
    trackBothTypes = pt.get<G4bool>("computation.trackBothTypes");
    nChargeToTrackTogether = pt.get<G4int>("computation.nChargeToTrackTogether");                //if more devide into subcharges
    initialDisplacement = pt.get<G4bool>("computation.useInitialDisplacement");
    chargeCloudSigma = 0;
    useDiffusionRepulsion = pt.get<G4bool>("computation.useDiffusionRepulsion");
    doTrapping = pt.get<G4bool>("computation.useTrapping");
    trappingTime = pt.get<G4double>("computation.trappingTime") * ns;

    // Default mobilities
    Default_Electron_Mobility = pt.get<G4double>(sensor + ".Default_Electron_Mobility") * cm2 / s;
    Default_Hole_Mobility     = pt.get<G4double>(sensor + ".Default_Hole_Mobility") * cm2 / s;
    Default_Electron_D        = pt.get<G4double>(sensor + ".Default_Electron_D") * cm2 / s;
    Default_Hole_D = pt.get<G4double>(sensor + ".Default_Hole_D") * cm2 / s;

    Default_LET = pt.get<G4double>(sensor + ".Default_LET") * um / keV;

    //mobility dependence on electric field
    Electron_AlphaField = pt.get<G4double>(sensor + ".Electron_AlphaField") * cm / s;
    Electron_ThetaField = pt.get<G4double>(sensor + ".Electron_ThetaField");
    Electron_TempNominal = pt.get<G4double>(sensor + ".Electron_TempNominal"); // [K]
//     Electron_Beta           = atof((configdata.operator[] ("Electron_Beta")).c_str());
    Electron_Beta = pt.get<G4double>(sensor + ".Electron_Beta");
    Electron_Saturation_Velocity = Electron_AlphaField * std::pow(1. + Electron_ThetaField * std::exp(Temperature / Electron_TempNominal), -1.0);
//     Electron_Saturation_Velocity = Electron_AlphaField*TMath::Power(1.0+Electron_ThetaField*TMath::Exp(Temperature/Electron_TempNominal),-1.0);

    Hole_AlphaField = pt.get<G4double>(sensor + ".Hole_AlphaField") * cm2 / s;
    Hole_ThetaField = pt.get<G4double>(sensor + ".Hole_ThetaField");
    Hole_TempNominal = pt.get<G4double>(sensor + ".Hole_TempNominal"); // [K]
    Hole_Beta = pt.get<G4double>(sensor + ".Hole_Beta");
    Hole_Saturation_Velocity = Hole_AlphaField * std::pow(1.0 + Hole_ThetaField * std::exp(Temperature / Hole_TempNominal), -1.0);

    //computational parameters
    pulsePrecision = pt.get<G4double>("computation.pulsePrecision") * ns;  //must be ns!
    maxPulseTime = pt.get<G4double>("computation.maxPulseTime") * ns;

    nPulseArrayElements = (G4int) maxPulseTime / pulsePrecision;
    ampResponseTime = pt.get<G4double>("computation.ampResponseTime") * ns;
    nAmpResponseElements = (G4int)(ampResponseTime / pulsePrecision) + nPulseArrayElements;
    cutDepthTracking = detectorThickness;
    writePeakToFile = pt.get<G4bool>("testing.writePeakToFile");
    useCSM = myDet->GetCsmMode();
    minEnergy = 0 * keV;

    DEBUG = pt.get<G4bool>("testing.debug");
    chargeCollectionMode = pt.get<G4bool>("testing.chargeCollectionMode");

    preamp = new PreampMedipix(pt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DigitizerWeightField::~DigitizerWeightField()
{
    delete[] potentialFine;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DigitizerWeightField::Digitize()
{
    // Get a pointer to the digiManager
    G4DigiManager *digMan = G4DigiManager::GetDMpointer();
    // Get the hit collection ID (defined in DetectorConstruction::ConstructSDandField())
    G4int MpxHitCollID = digMan->GetHitsCollectionID("SensorHitsCollection");
    // And fetch the Hits Collection
    const DetectorHitsCollection *hitCollection = static_cast<const DetectorHitsCollection *>(digMan->GetHitsCollection(MpxHitCollID));

    digitCollection = new MpxDigitCollection("DigitizerWeightField", "DigitsCollection");

    map<pair<G4int, G4int>, G4double * > *inducedPixelContent = new map<pair<G4int, G4int>, G4double * >;

    G4int nEntries = hitCollection->entries();
    if (nEntries == 0) return;

    G4int event = -1;


    //iterate over all entries in the list
    for (G4int itr = 0; itr < nEntries; itr++) {

        DetectorHit *hit = (*hitCollection)[itr];
        G4ThreeVector relPos;

        event = hit->GetEvent();  //FIXME if we run multiple particles per event later!!!

        //DebugData
//       relPos.set(55*um, 55*um, 50*um);
//       hit = new DetectorHit();
//       hit->Add(10* keV,0, relPos, 30, 30, 1, 1, 0*ns);
//       nEntries = 1;
//         G4cout << "col: " << hit->GetColumn() << " line: " << hit->GetLine() << " ener: " << hit->GetEdep()/keV << " pos: " << GetPixelCoordinates(hit->GetPosition()).x()/um << " " << GetPixelCoordinates(hit->GetPosition()).y()/um << " " << GetPixelCoordinates(hit->GetPosition()).z()/um << G4endl;
        //end Debug

        if (initialDisplacement == true) ComputeChargeCloudSigma(hit->GetEdep());
        G4double charge = FanoNoise((hit->GetEdep() / eV) / nElectronHolePairs); //(hit->GetEdep()/eV)/nElectronHolePairs;

        G4int nSubCharges = (G4int)charge / nChargeToTrackTogether;
        nSubCharges += 1; //resolve 0 case

        // devide charge in to subcharges for tracking
        for (G4int nQ = 0; nQ < nSubCharges; nQ++) {
            G4double subCharge = charge / nSubCharges;

            //position subcharge
            if (initialDisplacement == true) {
                relPos.setX(CLHEP::RandGauss::shoot(hit->GetPosition().x(), chargeCloudSigma));
                relPos.setY(CLHEP::RandGauss::shoot(hit->GetPosition().y(), chargeCloudSigma));
                relPos.setZ(CLHEP::RandGauss::shoot(hit->GetPosition().z(), chargeCloudSigma));
                if (relPos.z() > detectorThickness / 2) relPos.setZ(detectorThickness / 2);
                if (relPos.z() < -detectorThickness / 2) relPos.setZ(-detectorThickness / 2);
            } else {             
                relPos.set(hit->GetPosition().x(), hit->GetPosition().y(), hit->GetPosition().z());
            }

            //get pixel and relative position
            std::pair<G4int, G4int> tempPixel;
            tempPixel = GetPixelFromPosition(relPos);
            relPos = GetPixelCoordinates(relPos);

            //perform the drift
            ComputeDrift(relPos, tempPixel.first, tempPixel.second, hit->GetTime() / ns, subCharge, inducedPixelContent, typeToCollect, true);

            if (trackBothTypes) {
                //typeToCollect = !typeToCollect;
                ComputeDrift(relPos, tempPixel.first, tempPixel.second, hit->GetTime() / ns, subCharge, inducedPixelContent, !typeToCollect, false);
            }
        }
    }

    //write induced current peaks to file (before preamp)
    if (writePeakToFile == true) {
        std::ofstream myfile;
        map<pair<G4int, G4int>, G4double * >::iterator iPCItr;
        iPCItr = inducedPixelContent->begin();

        for (; iPCItr != inducedPixelContent->end() ; iPCItr++) {
            G4int k = (*iPCItr).first.first * 1000 + (*iPCItr).first.second;
            G4double *tempArray = (*iPCItr).second;

            char filename[64];
            std::sprintf(filename, "inducedCharge%d", k);

            myfile.open(filename);

            for (G4int it = 0; it < nPulseArrayElements; it++) {
                myfile << tempArray[it] << " ";
                //if((*iPCItr).first.first == 5 and (*iPCItr).first.second == 5) G4cout << "writeFile: " << tempArray[it] << G4endl;
            }
            myfile.close();
        }
    }


    // Electronics
    preamp->GetPixelResponse(inducedPixelContent, digitCollection, event);

    // Store the digit collection in the Digitizer Manager so it's available in the Run class
    StoreDigiCollection(digitCollection);

    //add pixel events to digitizer hit collection

    // TODO: The old export code runs from here, but has been disabled in favor of hdf5 export. This could (should?) be fixed
    //detector->AddPixelEvents(digitCollection);

    //print the contents of the digitCollection to the console
    /*
    G4int n = digitCollection->GetSize();
    for (G4int i = 0; i < n; i++) {
        G4cout << (*digitCollection)[i]->GetEvent() << " " << (*digitCollection)[i]->GetColumn() << " " << (*digitCollection)[i]->GetLine() << " " << (*digitCollection)[i]->GetEnergy() << G4endl;
    }
     */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DigitizerWeightField::CreatePotentialTable()
{
    G4String filename;
    
    if (detectorThickness == 200 * um) {
        if (pixelSize == 55 * um) {
            filename = "potential/3d_regular_grid_weighting_field55_200_fine.bin";
        }
        else{
            G4cout << "Error: no potential map for " << pixelSize / um << "um pixelsize, " << detectorThickness / um << "um sensor thickness" << G4endl;
            return;
        }
    } else if (detectorThickness == 300 * um) {
        if (pixelSize == 55 * um) {
            filename = "potential/3d_regular_grid_weighting_field55_300_fine.bin";
        } else if (pixelSize == 110 * um) {
            filename = "potential/3d_regular_grid_weighting_field110_300_fine.bin";
	} else if (pixelSize == 220 * um){
	    filename = "potential/3d_regular_grid_weighting_field220_300_fine.bin";
        } else {
            G4cout << "Error: no potential map for " << pixelSize / um << "um pixelsize, " << detectorThickness / um << "um sensor thickness" << G4endl;
            return;
        }

    } else if (detectorThickness == 500 * um) {
        if (pixelSize == 55 * um) {
            filename = "potential/3d_regular_grid_weighting_field55_500_fine.bin";
        } else if (pixelSize == 110 * um) {
            filename = "potential/3d_regular_grid_weighting_field110_500_fine.bin";
        } else {
            G4cout << "Error: no potential map for " << pixelSize / um << "um pixelsize, " << detectorThickness / um << "um sensor thickness" << G4endl;
            return;
        }
    } else if (detectorThickness == 1000 * um) {
        if (pixelSize == 55 * um) {
            filename = "potential/3d_regular_grid_weighting_field55_1000_fine.bin";
        } else if (pixelSize == 110 * um) {
            filename = "potential/3d_regular_grid_weighting_field110_1000_fine.bin";
        } else {
            G4cout << "Error: no potential map for " << pixelSize / um << "um pixelsize, " << detectorThickness / um << "um sensor thickness" << G4endl;
            return;
        }
    } else {
        G4cout << "Error: no potential map for " << pixelSize / um << "um pixelsize, " << detectorThickness / um << "um sensor thickness" << G4endl;
        return;
    }

    FILE *fp;
    G4double v;

    if ((fp = fopen(filename, "rb")) == NULL) {
        G4cout << "Error WFDigitizer. Could not open " << filename << G4endl;
        return;
    }

    //read grid parameters first
    fread((void *)(&v), sizeof(v), 1, fp);
    gridSizeFineXY  = v;
    fread((void *)(&v), sizeof(v), 1, fp);
    gridSizeFineZ   = v;
    fread((void *)(&v), sizeof(v), 1, fp);
    nGridPointsXY  = (G4int)v;
    fread((void *)(&v), sizeof(v), 1, fp);
    nGridPointsZ   = (G4int)v;

    potentialFine = new G4double[nGridPointsXY * nGridPointsXY * nGridPointsZ + nGridPointsXY * nGridPointsXY]();

    //fill the potential matrix
    for (G4int z = 1; z < nGridPointsZ + 1; z++) {
        for (G4int y = 0; y < nGridPointsXY; y++) {
            for (G4int x = 0; x < nGridPointsXY; x++) {
                fread((void *)(&v), sizeof(v), 1, fp);
                potentialFine[x + y * nGridPointsXY + z * nGridPointsXY * nGridPointsXY] = v;
            }
        }
    }

    //linear interpolation to get values behind the pixel
    for (G4int y = 0; y < nGridPointsXY; y++) {
        for (G4int x = 0; x < nGridPointsXY; x++) {
            potentialFine[x + y * nGridPointsXY] = potentialFine[x + y * nGridPointsXY + 1 * nGridPointsXY * nGridPointsZ]
                                                   + (potentialFine[x + y * nGridPointsXY + 1 * nGridPointsXY * nGridPointsZ]
                                                           - potentialFine[x + y * nGridPointsXY + 2 * nGridPointsXY * nGridPointsZ]);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DigitizerWeightField::GetInterpolatedPotential(G4ThreeVector position)
{
    //trilinear interpolation
    //conversion to um scale
    G4ThreeVector pos(position.x() / um, position.y() / um, position.z() / um);
    G4double c = 0;

    if(chargeCollectionMode == true){
        if (pos.z() == 0) {
            if (abs(pos.x()) <= (pixelSize / um) / 2 && abs(pos.y()) <= (pixelSize / um) / 2) c = 1; //on the pixel
            //}else c = 0;
        }
    } else {
        if (pos.z() == 0) {
            if (abs(pos.x()) <= (pixelSize / um) / 2 && abs(pos.y()) <= (pixelSize / um) / 2) c = 1; //on the pixel
            //}else c = 0;
        } else if (pos.z() >= cutDepthTracking / um) {
            c = 0;
        } else {
            //x
            G4double xmod = std::fmod(pos.x(), gridSizeFineXY);
            if (xmod < 0) xmod = gridSizeFineXY + xmod;
            G4double xd = xmod / gridSizeFineXY;
            G4int x0 = (int)(pos.x() / gridSizeFineXY);
            if (pos.x() < 0) {
                x0 += (nGridPointsXY / 2) - 1; //   49;
            } else {
                x0 += (nGridPointsXY / 2); //50;
            }

            //y
            G4double ymod = std::fmod(pos.y(), gridSizeFineXY);
            if (ymod < 0) ymod = gridSizeFineXY + ymod;
            G4double yd = ymod / gridSizeFineXY;
            G4int y0 = (int)(pos.y() / gridSizeFineXY);
            if (pos.y() < 0) {
                y0 += (nGridPointsXY / 2) - 1; //49;
            } else {
                y0 += (nGridPointsXY / 2); //50;
            }

            //z
            if (pos.z() < -(gridSizeFineZ / 2)) pos.setZ(0);
            G4double zmod = std::fmod(pos.z() - (gridSizeFineZ / 2), gridSizeFineZ);
            G4double zd = zmod / gridSizeFineZ;
            G4int z0 = (int)(((pos.z() - gridSizeFineZ / 2) / gridSizeFineZ));

            z0 = z0 + 1;

            if (x0 < 0 || x0+1 > nGridPointsXY || y0 < 0 || y0+1 > nGridPointsXY || z0+1 > nGridPointsZ) {
                //G4cout << "Error in GetInterpolatedPotential: Index out of range. " << " x " << x0 << " y " << y0 << "z " << z0 << G4endl;
                c = 0;
            } else {

                //conversion to interpolated matrix with values behind detector

                G4double c00 = GetPotential(x0, y0, z0) * (1 - xd) + GetPotential(x0 + 1, y0, z0)    * xd;
                G4double c10 = GetPotential(x0, y0 + 1, z0) * (1 - xd) + GetPotential(x0 + 1, y0 + 1, z0)  * xd;
                G4double c01 = GetPotential(x0, y0, z0 + 1) * (1 - xd) + GetPotential(x0 + 1, y0, z0 + 1)  * xd;
                G4double c11 = GetPotential(x0, y0 + 1, z0 + 1) * (1 - xd) + GetPotential(x0 + 1, y0 + 1, z0 + 1) * xd;

                G4double c0 = c00 * (1 - yd) + c10 * yd;
                G4double c1 = c01 * (1 - yd) + c11 * yd;
                c = c0 * (1 - zd) + c1 * zd;
            }
        }
    }
    return c;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DigitizerWeightField::GetPotential(G4int x, G4int y, G4int z)
{
    G4double tempPot = potentialFine[x + y * nGridPointsXY + z * nGridPointsXY * nGridPointsXY];
    return tempPot;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DigitizerWeightField::ComputeElectricField(G4double /*x*/, G4double /*y*/, G4double z)
{
    electricFieldX = 0;
    electricFieldY = 0;
    if (material == "G4_Si") {
        ComputeElectricField1D(z);
    }
    else if (material == "G4_CADMIUM_TELLURIDE")
    {
      //electricFieldZ;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DigitizerWeightField::ComputeElectricField1D(G4double z)
{
    //electric field for overdepleted sensor
    if (z > depletedDepth) {
        electricFieldZ = 0;
    } else {
        electricFieldZ = 2 * depletionVoltage / detectorThickness * (1 - z/detectorThickness) + (biasVoltage - depletionVoltage) / detectorThickness;
        if(electricFieldZ < 0) electricFieldZ = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DigitizerWeightField::GetElectricFieldNorm(G4double /*x*/, G4double /*y*/, G4double /*z*/)
{
    return std::sqrt(electricFieldX * electricFieldX + electricFieldY * electricFieldY + electricFieldZ * electricFieldZ);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4double>  DigitizerWeightField::ComputeDrift(G4ThreeVector position, G4int column, G4int line, G4double time, G4double charge, map<pair<G4int, G4int>, G4double * > *pInducedPixelContent, G4bool particle, G4bool primType)
{
    G4double driftTime = 0;
    G4double globalTime = time;
    G4double dt = pulsePrecision - fmod(time, pulsePrecision); // set first dt to reach grid, then gridsize
    G4double xtemp = position.x();
    G4double ytemp = position.y();
    G4double ztemp = position.z();
    vector<G4double> step(4, 0);
    int counter = 0;
    G4ThreeVector lastPosition(xtemp, ytemp, ztemp);
    G4double potLastPosition[3][3] = {{0}};
    G4double potNewPosition[3][3] = {{0}};
    G4ThreeVector newPosition(0., 0., 0.);
    G4double *pulseArray[3][3];
    G4double chargeAfterTrapping = charge;
    G4double inducedCharge = 0;
    G4double sigma = std::max(1 * um, chargeCloudSigma); //initial sigma but at least 1um

    //initialize potential
    for (G4int i = -1; i <= 1; i++) {
        for (G4int j = -1; j <= 1; j++) {
            G4ThreeVector tempPosition(xtemp - (pixelSize * i), ytemp - (pixelSize * j), ztemp);
            potLastPosition[i + 1][j + 1] =  GetInterpolatedPotential(tempPosition);
        }
    }

    //data structure
    std::pair<G4int, G4int> extrapixel;
    extrapixel.first = column;
    extrapixel.second = line;
    std::pair<G4int, G4int> tempPixel;
    map<pair<G4int, G4int>, G4double * >::iterator pCItr;

    //search for pixel or create
    for (G4int i = -1; i <= 1; i++) {
        for (G4int j = -1; j <= 1; j++) {
            tempPixel = extrapixel;
            tempPixel.first     += i;
            tempPixel.second    += j;

            pCItr = pInducedPixelContent->find(tempPixel);
            if (pCItr == pInducedPixelContent->end()) {
                //check if outside of detector
                if (tempPixel.first >= 0 && tempPixel.first < nPixX && tempPixel.second >= 0 && tempPixel.second < nPixY) {
                    pulseArray[i + 1][j + 1] = new G4double[nPulseArrayElements];
                    std::fill_n(pulseArray[i + 1][j + 1], nPulseArrayElements, 0);
                    pInducedPixelContent->operator[](tempPixel) = pulseArray[i + 1][j + 1];
                } else {
                    pulseArray[i + 1][j + 1] = NULL;
                }
            } else {
                pulseArray[i + 1][j + 1] = (*pCItr).second ;
            }
        }
    }

    //perform the stepping
    while (ztemp > 0 && ztemp < cutDepthTracking) {

        driftTime += dt;
        globalTime += dt;
        if (globalTime >= maxPulseTime) break;
        if (doTrapping == true && typeToCollect == 1) chargeAfterTrapping = charge * std::exp(-driftTime / trappingTime);

        step = RKF45Integration(xtemp, ytemp, ztemp, dt, particle);
        sigma = ComputeDiffusionRepulsionSigma(dt, sigma, chargeAfterTrapping, particle);
        if (primType) {
            xtemp += step[0];
            ytemp += step[1];
            ztemp += step[2];
            //G4cout << "primary type" << ztemp/um << G4endl;
        } else { //track away from pixel
            xtemp -= step[0];
            ytemp -= step[1];
            ztemp -= step[2];
            //G4cout << "secondary type" << ztemp/um << G4endl;
        }
        if (useDiffusionRepulsion == true) {
            xtemp = CLHEP::RandGauss::shoot(xtemp, sigma);
            ytemp = CLHEP::RandGauss::shoot(ytemp, sigma);
            ztemp = CLHEP::RandGauss::shoot(ztemp, sigma);
        }
        if (ztemp > detectorThickness) ztemp = detectorThickness;
        if (ztemp < 0) ztemp = 0;
        //if (GetElectricFieldNorm(xtemp, ytemp, ztemp) == 0)ztemp = detectorThickness;
        newPosition.set(xtemp, ytemp, ztemp);

        //additional timestep to make sure the float is casted correctly
        G4int index = (globalTime + 0.0001) / pulsePrecision;

        for (G4int i = -1; i <= 1; i++) {
            for (G4int j = -1; j <= 1; j++) {
                if (pulseArray[i + 1][j + 1] != NULL) {
                    G4ThreeVector tempPosition(xtemp - (pixelSize * i), ytemp - (pixelSize * j), ztemp);
                    potNewPosition[i + 1][j + 1] =  GetInterpolatedPotential(tempPosition);

                    inducedCharge = chargeAfterTrapping * (potNewPosition[i + 1][j + 1] - potLastPosition[i + 1][j + 1]); //Ramo-Shockley

                    //if(i == 0 and j == 0) G4cout << "induced Charge: " << inducedCharge << " index: " << index << G4endl;

                    //store charge
                    if (primType) {
                        pulseArray[i + 1][j + 1][index] += inducedCharge;
                    } else {
                        pulseArray[i + 1][j + 1][index] -= inducedCharge;
                    }
                }
            }
        }

        //store potential in lastPosition
        lastPosition = newPosition;
        for (G4int i = 0; i <= 3; i++) {
            for (G4int j = 0; j <= 3; j++) {
                potLastPosition[i][j] = potNewPosition[i][j];
            }
        }

        dt = pulsePrecision;
        counter++;
        //G4cout << "drift Time: " << driftTime << " iterations: " << counter << " ztemp: " << ztemp / um << " x: " << xtemp/um << " y: " << ytemp/um << " pot: " << potLastPosition[1][1] << " efield: " << electricFieldZ << " particle " << particle <<G4endl;
    }

    std::vector<G4double> output(4);
    output[0] = xtemp;
    output[1] = ytemp;
    output[2] = ztemp;
    output[3] = driftTime;
    //G4cout << "iterations: " << counter << G4endl;
    return output;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4double>  DigitizerWeightField::RKF45Integration(G4double x, G4double y, G4double z, G4double dt, G4bool particle)
{
// This function transport using Euler integration, for field (Ex,Ey,Ez),
// considered constant over time dt. The movement equation are those
// of charges in semi-conductors, sx= mu*E*dt;;
    double k1x, k2x, k3x, k4x, k5x, k6x;
    double k1y, k2y, k3y, k4y, k5y, k6y;
    double k1z, k2z, k3z, k4z, k5z, k6z;
    double dx, dy, dz;

    ComputeElectricField(x, y, z);

    k1x = -MobilityElectronHole(x, y, z, particle) * electricFieldX * dt;
    k1y = -MobilityElectronHole(x, y, z, particle) * electricFieldY * dt;
    k1z = -MobilityElectronHole(x, y, z, particle) * electricFieldZ * dt;

    ComputeElectricField(x + k1x / 4, y + k1y / 4, z + k1z / 4);

    k2x = -1 * MobilityElectronHole(x + k1x / 4, y + k1y / 4, z + k1z / 4, particle) * electricFieldX * dt;
    k2y = -MobilityElectronHole(x + k1x / 4, y + k1y / 4, z + k1z / 4, particle) * electricFieldY * dt;
    k2z = -MobilityElectronHole(x + k1x / 4, y + k1y / 4, z + k1z / 4, particle) * electricFieldZ * dt;

    ComputeElectricField(x + (9. / 32)*k2x + (3. / 32)*k1x, y + (9. / 32)*k2y + (3. / 32)*k1y, z + (9. / 32)*k2z + (3. / 32)*k1z);

    k3x = -MobilityElectronHole(x + (9. / 32) * k2x + (3. / 32) * k1x, y + (9. / 32) * k2y + (3. / 32) * k1y, z + (9. / 32) * k2z + (3. / 32) * k1z, particle) * electricFieldX * dt;
    k3y = -MobilityElectronHole(x + (9. / 32) * k2x + (3. / 32) * k1x, y + (9. / 32) * k2y + (3. / 32) * k1y, z + (9. / 32) * k2z + (3. / 32) * k1z, particle) * electricFieldY * dt;
    k3z = -MobilityElectronHole(x + (9. / 32) * k2x + (3. / 32) * k1x, y + (9. / 32) * k2y + (3. / 32) * k1y, z + (9. / 32) * k2z + (3. / 32) * k1z, particle) * electricFieldZ * dt;

    ComputeElectricField(x - (7200. / 2197)*k2x + (1932. / 2197)*k1x + (7296. / 2197)*k3x, y - (7200. / 2197)*k2y + (1932. / 2197)*k1y + (7296. / 2197)*k3y, z - (7200. / 2197)*k2z + (1932. / 2197)*k1z + (7296. / 2197)*k3z);

    k4x = -MobilityElectronHole(x - (7200. / 2197) * k2x + (1932. / 2197) * k1x + (7296. / 2197) * k3x, y - (7200. / 2197) * k2y + (1932. / 2197) * k1y + (7296. / 2197) * k3y, z - (7200. / 2197) * k2z + (1932. / 2197) * k1z + (7296. / 2197) * k3z, particle) * electricFieldX * dt;
    k4y = -MobilityElectronHole(x - (7200. / 2197) * k2x + (1932. / 2197) * k1x + (7296. / 2197) * k3x, y - (7200. / 2197) * k2y + (1932. / 2197) * k1y + (7296. / 2197) * k3y, z - (7200. / 2197) * k2z + (1932. / 2197) * k1z + (7296. / 2197) * k3z, particle) * electricFieldY * dt;
    k4z = -MobilityElectronHole(x - (7200. / 2197) * k2x + (1932. / 2197) * k1x + (7296. / 2197) * k3x, y - (7200. / 2197) * k2y + (1932. / 2197) * k1y + (7296. / 2197) * k3y, z - (7200. / 2197) * k2z + (1932. / 2197) * k1z + (7296. / 2197) * k3z, particle) * electricFieldZ * dt;

    ComputeElectricField(x - (8)*k2x + (439. / 216)*k1x + (3680. / 513)*k3x - (845. / 4104)*k4x, y - (8)*k2y + (439. / 216)*k1y + (3680. / 513)*k3y - (845. / 4104)*k4y, z - (8)*k2z + (439. / 216)*k1z + (3680. / 513)*k3z - (845. / 4104)*k4z);

    k5x = -MobilityElectronHole(x - (8) * k2x + (439. / 216) * k1x + (3680. / 513) * k3x - (845. / 4104) * k4x, y - (8) * k2y + (439. / 216) * k1y + (3680. / 513) * k3y - (845. / 4104) * k4y, z - (8) * k2z + (439. / 216) * k1z + (3680. / 513) * k3z - (845. / 4104) * k4z, particle) * electricFieldX * dt;
    k5y = -MobilityElectronHole(x - (8) * k2x + (439. / 216) * k1x + (3680. / 513) * k3x - (845. / 4104) * k4x, y - (8) * k2y + (439. / 216) * k1y + (3680. / 513) * k3y - (845. / 4104) * k4y, z - (8) * k2z + (439. / 216) * k1z + (3680. / 513) * k3z - (845. / 4104) * k4z, particle) * electricFieldY * dt;
    k5z = -MobilityElectronHole(x - (8) * k2x + (439. / 216) * k1x + (3680. / 513) * k3x - (845. / 4104) * k4x, y - (8) * k2y + (439. / 216) * k1y + (3680. / 513) * k3y - (845. / 4104) * k4y, z - (8) * k2z + (439. / 216) * k1z + (3680. / 513) * k3z - (845. / 4104) * k4z, particle) * electricFieldZ * dt;

    ComputeElectricField(x + (2)*k2x - (8. / 27)*k1x - (3544. / 2565)*k3x - (1859. / 4104)*k4x - (11. / 40)*k5x,
                         y + (2)*k2y - (8. / 27)*k1y - (3544. / 2565)*k3y - (1859. / 4104)*k4y - (11. / 40)*k5y,
                         z + (2)*k2z - (8. / 27)*k1z - (3544. / 2565)*k3z - (1859. / 4104)*k4z - (11. / 40)*k5z);

    k6x = -MobilityElectronHole(x + (2) * k2x - (8. / 27) * k1x - (3544. / 2565) * k3x - (1859. / 4104) * k4x - (11. / 40) * k5x,
                                y + (2) * k2y - (8. / 27) * k1y - (3544. / 2565) * k3y - (1859. / 4104) * k4y - (11. / 40) * k5y,
                                z + (2) * k2z - (8. / 27) * k1z - (3544. / 2565) * k3z - (1859. / 4104) * k4z - (11. / 40) * k5z, particle) * electricFieldX * dt;
    k6y = -MobilityElectronHole(x + (2) * k2x - (8. / 27) * k1x - (3544. / 2565) * k3x - (1859. / 4104) * k4x - (11. / 40) * k5x,
                                y + (2) * k2y - (8. / 27) * k1y - (3544. / 2565) * k3y - (1859. / 4104) * k4y - (11. / 40) * k5y,
                                z + (2) * k2z - (8. / 27) * k1z - (3544. / 2565) * k3z - (1859. / 4104) * k4z - (11. / 40) * k5z, particle) * electricFieldY * dt;
    k6z = -MobilityElectronHole(x + (2) * k2x - (8. / 27) * k1x - (3544. / 2565) * k3x - (1859. / 4104) * k4x - (11. / 40) * k5x,
                                y + (2) * k2y - (8. / 27) * k1y - (3544. / 2565) * k3y - (1859. / 4104) * k4y - (11. / 40) * k5y,
                                z + (2) * k2z - (8. / 27) * k1z - (3544. / 2565) * k3z - (1859. / 4104) * k4z - (11. / 40) * k5z, particle) * electricFieldZ * dt;

    dx = ((16. / 135) * k1x + (6656. / 12825) * k3x + (28561. / 56430) * k4x - (9. / 50) * k5x + (2. / 55) * k6x);
    dy = ((16. / 135) * k1y + (6656. / 12825) * k3y + (28561. / 56430) * k4y - (9. / 50) * k5y + (2. / 55) * k6y);
    dz = ((16. / 135) * k1z + (6656. / 12825) * k3z + (28561. / 56430) * k4z - (9. / 50) * k5z + (2. / 55) * k6z);

    std::vector<G4double> newpoint(4, 0);
    newpoint[0] = (dx);
    newpoint[1] = (dy);
    newpoint[2] = (dz);
    newpoint[3] = 0; //no error handling

    return newpoint;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double  DigitizerWeightField::MobilityElectronHole(G4double /*x*/, G4double /*y*/, G4double /*z*/, G4bool particle)
{
    G4double parElectricField = electricFieldZ;//GetElectricFieldNorm(x,y,z);

    if (particle == 0) {

        mobility = Default_Electron_Mobility *
                   std::pow((1.0 / (1. + std::pow((Default_Electron_Mobility * parElectricField) /
                                    Electron_Saturation_Velocity, Electron_Beta))), 1.0 / Electron_Beta);

    } else if (particle == 1) {
        mobility = Default_Hole_Mobility *
                   std::pow((1.0 / (1. + std::pow((Default_Hole_Mobility * parElectricField) /
                                    Hole_Saturation_Velocity, Hole_Beta))), 1.0 / Hole_Beta);

        //mobility = - mobility;  //relative minus for orientation
    }
    return mobility;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector DigitizerWeightField::GetPixelCoordinates(G4ThreeVector pos)
{
    //TODO should pixelSize be integer in DetectorConstruction?
    G4double x = pos.x();
    G4double y = pos.y();

    G4ThreeVector pixelPos(0., 0., 0.);

    //even number pixels
    if (nPixX % 2 == 0 && nPixY % 2 == 0) {
        x = std::fmod(x, pixelSize);
        y = std::fmod(y, pixelSize);
        if (x >= 0) {
            x = x - pixelSize / 2.;
        } else {
            x = x + pixelSize / 2.;
        }

        if (y >= 0) {
            y = y - pixelSize / 2.;
        } else {
            y = y + pixelSize / 2.;
        }
        //relative minus in z caused by inverse z direction
        pixelPos.set(x, y, detectorThickness / 2. - pos.z());

        //uneven number pixels
    } else {
        x = std::fmod(x, pixelSize);
        y = std::fmod(y, pixelSize);
        if (x >= 0) {
            if (x >= pixelSize / 2) x = x - pixelSize;
        } else {
            if (x < -pixelSize / 2) x = x + pixelSize;
        }

        if (y >= 0) {
            if (y >= pixelSize / 2) y = y - pixelSize;
        } else {
            if (y < -pixelSize / 2) y = y + pixelSize;
        }
        //relative minus in z caused by inverse z direction
        pixelPos.set(x, y, detectorThickness / 2. - pos.z());
    }
    return pixelPos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::pair<G4int, G4int> DigitizerWeightField::GetPixelFromPosition(G4ThreeVector pos)
{
    std::pair<G4int, G4int> pixel;
    if (nPixX % 2 == 0 || nPixY % 2 == 0) {
        pixel.first = (G4int)((pos.x() + nPixX / 2 * pixelSize) / pixelSize);
        pixel.second = (G4int)((pos.y() + nPixY / 2 * pixelSize) / pixelSize);
        //uneven number pixel
    } else {
        pixel.first = (G4int)((pos.x() + ((nPixX - 1) / 2 + 0.5) * pixelSize) / pixelSize) + 1;
        pixel.second = (G4int)((pos.y() + ((nPixY - 1) / 2 + 0.5) * pixelSize) / pixelSize) + 1;
    }
    return pixel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DigitizerWeightField::FanoNoise(G4double charge)
{
    G4double sigma = std::sqrt(charge * fanoFactor);
    return CLHEP::RandGauss::shoot(charge, sigma);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DigitizerWeightField::ComputeChargeCloudSigma(G4double energy)
{
    chargeCloudSigma = Default_LET * energy * (1 - (0.98 / (1 + 0.003 * (1 / keV) * energy)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double DigitizerWeightField::ComputeDiffusionRepulsionSigma(G4double time, G4double sigma, G4double charge, G4bool particle)
{
    G4double d;
    G4double defaultMobility;
    switch (particle) {
    case 0:
        d = Default_Electron_D;
        defaultMobility = Default_Electron_Mobility;
        break;
    case 1:
        d = Default_Hole_D;
        defaultMobility = Default_Hole_Mobility;
        break;
    }

    //return sqrt(2 * d * time);
    return sqrt(2.*(d + (defaultMobility * std::abs(charge)) / (sigma * Default_Relative_Permittivity) * 1.354 * std::pow(10, -10) * m) * time);
}

