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
/// \file src/MpxDetector.cc
/// \brief Implementation of the MpxDetector class
/// Collects information from different detector typse (Medipix3RX, Timepix
/// and writes results from digitizer to either sparse or ascii matrix
///

#include "G4AutoLock.hh"
#include "MpxDetector.hh"
#include "DetectorConstructionBase.hh"
//#include "config.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include "PrimaryGeneratorMessenger.hh"

#include <iostream>
#include <fstream>
#include <ctime>
#include "math.h"
#include <boost/concept_check.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <ExportMgr.hh>
#include <RunAction.hh>

// #include "../include/MpxDetector.hh"

#include "CLHEP/Random/RandGauss.h"

MpxDetector *MpxDetector::instance = 0;

MpxDetector::MpxDetector():
    writeCounter(0),
    nPixel(0)
{
    G4RunManager *fRM = G4RunManager::GetRunManager();
    DetectorConstructionBase *myDet = (DetectorConstructionBase *)(fRM->GetUserDetectorConstruction());
    nPixel = myDet->GetNbPixels();


    detectorType = myDet->GetDetectorType();
    // 0 = Medipix3RX
    // 1 = Timepix
    // 2 = Dosepix
    csmMode = myDet->GetCsmMode();


    //Medipix3RX
    MpxThreshold1 = 20.0 * keV;
    MpxThreshold2 = 20.0 * keV;
    MpxThreshold3 = 20.0 * keV;
    MpxThreshold4 = 20.0 * keV;

    //pixel matrix
    MpxPixelMat1 = new G4int[nPixel * nPixel]();
    MpxPixelMat2 = new G4int[nPixel * nPixel]();
    MpxPixelMat3 = new G4int[nPixel * nPixel]();
    MpxPixelMat4 = new G4int[nPixel * nPixel]();


    //Timepix configuration
    tpxMode = myDet->GetTpxMode();

    // emty frame with nxn pixels
    TpxPixelMat = new G4double[nPixel * nPixel]();
    TpxThlDisp = new G4double[nPixel * nPixel]();

    //Output
    frameCounter = 0;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
MpxDetector::~MpxDetector()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace
{
G4Mutex singletonMutex = G4MUTEX_INITIALIZER;
}

MpxDetector *MpxDetector::GetInstance()
{
    G4AutoLock l(&singletonMutex);
    if (!instance) {
        instance = new MpxDetector();
    }
    return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace
{
G4Mutex AddDataMutex = G4MUTEX_INITIALIZER;
}

void MpxDetector::AddPixelEvents(MpxDigitCollection *collection)
{
    G4AutoLock l2(&AddDataMutex);

    for (G4int i = 0; i < (int)collection->GetSize(); i++) {
        Digit *digit = (*collection)[i];
        G4int pixelIndex = (digit->GetColumn()) * nPixel + (digit->GetLine());

        // set the 4 threshold in the Medipix3RX detector
        if (detectorType == 0) {
            //MPX3RX
            if (digit->GetEnergy() >= MpxThreshold1 / keV) MpxPixelMat1[pixelIndex] += 1;
            if (digit->GetEnergy() >= MpxThreshold2 / keV) MpxPixelMat2[pixelIndex] += 1;
            if (digit->GetEnergy() >= MpxThreshold3 / keV) MpxPixelMat3[pixelIndex] += 1;
            if (digit->GetEnergy() >= MpxThreshold4 / keV) MpxPixelMat4[pixelIndex] += 1;

            //Sparse stuff
            struct snglEvent newSnglEvent =  {(uint32_t)digit->GetEvent(), (uint32_t)digit->GetColumn(),
                (uint32_t)digit->GetLine(), (G4double)digit->GetEnergy(), (G4double)digit->GetToT(), (G4double)digit->GetToA()};
            sparseList.push_back(newSnglEvent);
            digit->~Digit();
        }

        //Timepix
        //Timepix operation modes 0 = PC, 1 = ToT, 2 = ToA
        if (detectorType == 1 || detectorType == 3) {

            if (tpxMode == 0) {
                //G4cout << "Counting mode" << G4endl;
                //G4cout << "DEBUG: TpxThreshold is " << TpxThreshold /keV << G4endl;
                //G4cout << "DEBUG: Photon energy is " << digit->GetEnergy()  << G4endl;
                //if (digit->GetEnergy() >= (TpxThreshold / keV)) {
                    TpxPixelMat[pixelIndex] += 1;
                    //G4cout << "Counting" << G4endl;
                //}
            } else if (tpxMode == 1) {
                //G4cout << "ToT mode" << G4endl;
                //if (digit->GetEnergy() >= (TpxThreshold / keV)) {
                    TpxPixelMat[pixelIndex] = digit->GetEnergy();
                //}
            }

            //Sparse stuff
            //if (digit->GetEnergy() >= (TpxThreshold / keV)) {
            struct snglEvent newSnglEvent =  {(uint32_t)digit->GetEvent(), (uint32_t)digit->GetColumn(),
                        (uint32_t)digit->GetLine(), (G4double)digit->GetEnergy(), (G4double)digit->GetToT(), (G4double)digit->GetToA()};
            sparseList.push_back(newSnglEvent);
	       //G4cout << "MpxDetector thdisp in pixel: " << TpxThlDisp[pixelIndex] << G4endl;
           //}
           digit->~Digit();
        } //end Timpepix


        //Dosepix
        if (detectorType == 2) {
            //Sparse stuff
            //Mybe implement the binning later

            G4double  eDep = digit->GetEnergy();

            //Transform to ToT hardcoded values
            G4double a = 2.66;
            G4double b = 126.13;
            G4double c = 375.02;
            G4double t = 0.5764;

            G4double tot = a * eDep + b - c / (eDep - t);

            if ((eDep > 2) && (tot > 0)) {
                tot = ceil(tot);
                //back to energy for storing in sparse file
                //G4double newEdep = (t * a - b + tot + sqrt(pow((b + t * a - tot), 2) + 4 * a * c)) / (2 * a);

                struct snglEvent newSnglEvent =  {(uint32_t)digit->GetEvent(), (uint32_t)digit->GetColumn(), (uint32_t)digit->GetLine(), (G4double)eDep, (G4double)tot, digit->GetToA()};
                sparseList.push_back(newSnglEvent);
                digit->~Digit();

            }
        } // end Dosepix
    }
    writeCounter++;
    if (writeCounter == PIXELS_CHUNK_SIZE) {
        WriteSparse();
        writeCounter = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace
{
G4Mutex WriteMutex = G4MUTEX_INITIALIZER;
}

void MpxDetector::WriteSparse()
{
    G4AutoLock l3(&WriteMutex);

    G4cout << "Writing sparse pixels output per event. Number of digits: " << sparseList.size() << G4endl;

    // Access to detector parameters
    G4RunManager *fRM = G4RunManager::GetRunManager();
    DetectorConstructionBase *myDet = (DetectorConstructionBase *)(fRM->GetUserDetectorConstruction());
    G4String sfname = myDet->GetSparseOutputFilename();

    //binary sparse export
    if (sfname != "") {
        std::ofstream myFile(sfname, std::ios::out | std::ios::binary | std::ios::app);
        for (std::list<snglEvent>::const_iterator iterator = sparseList.begin(); iterator != sparseList.end(); ++iterator) {
            struct snglEvent newSnglEvent = *iterator;
            //G4cout << newSnglEvent.col << "------" << sizeof(newSnglEvent) << G4endl;
            //G4cout << sizeof(newSnglEvent.event) << sizeof(newSnglEvent.col) << sizeof(newSnglEvent.line) << sizeof(newSnglEvent.energy) << sizeof(newSnglEvent) <<G4endl;
            myFile.write((char *)&newSnglEvent.event, 4);
            myFile.write((char *)&newSnglEvent.col, 4);
            myFile.write((char *)&newSnglEvent.line, 4);
            myFile.write((char *)&newSnglEvent.energy, 8);
            myFile.write((char *)&newSnglEvent.tot, 8);
            myFile.write((char *)&newSnglEvent.toa, 8);
        }

        myFile.close();
    }

#ifdef WITH_HDF5
    //RunAction * uSR = (RunAction * ) fRM->GetUserRunAction();
    //uSR->getExportManager()->WritePixels(sparseList);
#endif
    sparseList.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MpxDetector::WriteFrame()
{
    // Write the frame files. Should only be done at the end of each run.

    G4AutoLock l3(&WriteMutex);


    G4cout << "frameCounter: " << frameCounter << G4endl;

    //Acces to detector parameters and ascii file name
    G4RunManager *fRM = G4RunManager::GetRunManager();
    DetectorConstructionBase *myDet = (DetectorConstructionBase *)(fRM->GetUserDetectorConstruction());
    G4String fname = myDet->GetOutputFilename();

    //Frame output
    if (fname != "" ){
	//Medipix3RX: write frame matrix to ascii file
	if (detectorType == 0) {
	    //ASCII pixel matrix
	    std::ofstream myfile;

	    myfile.open(fname);

	    for (G4int i = 0; i < nPixel; i++) {
		for (G4int j = 0; j < nPixel; j++) {
		    myfile << MpxPixelMat1[i * nPixel + j] << " ";
		    myfile << MpxPixelMat2[i * nPixel + j] << " ";
		    if (j == nPixel - 1) myfile << "\n";
		}
		for (G4int j = 0; j < nPixel; j++) {
		    myfile << MpxPixelMat3[i * nPixel + j] << " ";
		    myfile << MpxPixelMat4[i * nPixel + j] << " ";
		    if (j == nPixel - 1) myfile << "\n";
		}
	    }
	    myfile.close();
	    G4cout << "MpxDetector (Medipix) writing file" << G4endl;
	}


	//Timepix: write frame matrix to ascii file
	if (detectorType == 1) {
	  
	  //Binary Matrix single File
	  G4bool sparse = false;
	  if (sparse == false){
	    std::ios_base::openmode flag;
	    if(frameCounter == 0){
	      flag = (std::ios::out | std::ios::binary | std::ios::trunc);
	    }else{
	      flag = (std::ios::out | std::ios::binary | std::ios::app);
	    }
	    std::ofstream myFile(fname, flag);
	    myFile.write((char *)TpxPixelMat, nPixel*nPixel*sizeof(G4double));
	    
	    G4cout << "MpxDetector (Timepix) writing file" << G4endl;
	  }//if sparse == false
	  else {
	    // Sparse frame
	  }
	}
    frameCounter++;
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MpxDetector::ResetPixelMatrix()
{
    //FIXME Check for detector type!!!
    for (G4int i = 0; i < nPixel * nPixel; i++)
        TpxPixelMat[i] = 0;

    for (G4int i = 0; i < nPixel * nPixel; i++) {
        MpxPixelMat1[i] = 0;
        MpxPixelMat2[i] = 0;
        MpxPixelMat3[i] = 0;
        MpxPixelMat4[i] = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MpxDetector::ReInitMatrix()
{
    //Needed to have the right size of the pixel matrix

    G4RunManager *fRM = G4RunManager::GetRunManager();
    DetectorConstructionBase *myDet = (DetectorConstructionBase *)(fRM->GetUserDetectorConstruction());
    nPixel = myDet->GetNbPixels();

    delete[] MpxPixelMat1;
    delete[] MpxPixelMat2;
    delete[] MpxPixelMat3;
    delete[] MpxPixelMat4;

    //pixel matrix
    MpxPixelMat1 = new G4int[nPixel * nPixel]();
    MpxPixelMat2 = new G4int[nPixel * nPixel]();
    MpxPixelMat3 = new G4int[nPixel * nPixel]();
    MpxPixelMat4 = new G4int[nPixel * nPixel]();


    delete[] TpxPixelMat;
    TpxPixelMat = new G4double[nPixel * nPixel]();

    delete[] TpxThlDisp;
    TpxThlDisp = new G4double[nPixel * nPixel]();

    SetThresholdDispersion();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MpxDetector::SetThresholdDispersion()
{
    //load .ini file with configuration data
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("DetectorConfig.ini", pt);

    //sensor properties
    G4RunManager *fRM = G4RunManager::GetRunManager();
    DetectorConstructionBase* myDet = (DetectorConstructionBase *)(fRM->GetUserDetectorConstruction());
    G4String material =  myDet->GetSensorMaterial()->GetName();
    G4String sensor;
    if (material == "G4_Si")
        sensor = "sensor_silicon";
    else if (material == "G4_CADMIUM_TELLURIDE")
    {
        sensor = "sensor_cdte";
    }

    G4double disp = 0;
    G4double nElectronHolePairs = pt.get<G4double>(sensor + ".nElectronHolePairs");

    if (detectorType == 0){
        if (csmMode == true){
            disp = pt.get<G4double>("chip_medipix3rx.csm_threshold_dispersion") * nElectronHolePairs * eV;
        } else {
            disp = pt.get<G4double>("chip_medipix3rx.spm_threshold_dispersion") * nElectronHolePairs * eV;
        }
    }else if (detectorType == 1){
        disp = pt.get<G4double>("chip_timepix.threshold_dispersion") * nElectronHolePairs * eV;
    }else if(detectorType == 2){
        disp = pt.get<G4double>("chip_dosepix.threshold_dispersion") * nElectronHolePairs * eV;
    }

    for (G4int i = 0; i < nPixel * nPixel; i++) {
        //G4cout << "MpxDetector: Setting Threshold Dispersion " << disp << " i: " << i << " nPixel: " << nPixel*nPixel << G4endl;
        TpxThlDisp[i] = CLHEP::RandGauss::shoot(0, disp/keV); //sigma of thl dispersion
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MpxDetector::WriteSimulationSettings()
{
    //Acces to detector parameters and ascii file name
    G4RunManager *fRM = G4RunManager::GetRunManager();
    DetectorConstructionBase *myDet = (DetectorConstructionBase *)(fRM->GetUserDetectorConstruction());

    // get current date and time
    std::time_t rawtime;
    std::tm *timeinfo;
    char buffer [80];

    std::time(&rawtime);
    timeinfo = std::localtime(&rawtime);

    std::strftime(buffer, 80, "%Y-%m-%d-%H-%M-%S", timeinfo);

    G4String date(buffer);
    G4String fname = myDet->GetConfigFilename();
    
    //No config filename specified use date and time    
    if (fname == "")
    {
      fname = date + "run_config.ini";
    }
    else
    {
      fname = date + "_config.ini";
    }
    using boost::property_tree::ptree;

    ptree pt;
    //write configuration parameters
    pt.put("sensor.nbPixels", myDet->GetNbPixels());
    pt.add("sensor.PixelSize", myDet->GetPixelSize());
    pt.add("sensor.SensorThickness", myDet->GetSensorThickness());
    pt.add("sensor.SensorMaterial", myDet->GetSensorMaterial()->GetName());
    pt.add("digitizer.DigitizerName", myDet->GetDigitizerName());
    G4String sdetectorType;
    if (detectorType == 0) {
      sdetectorType = "Medipix3RX";
    } else if (detectorType == 1) {
      sdetectorType = "Timepix";
    } else {
      sdetectorType = "";
    }
    pt.add("detector.DetectorType", sdetectorType);

    pt.add("detector.tpxMode", myDet->GetTpxMode());

    PrimaryGeneratorMessenger *prmMessenger = PrimaryGeneratorMessenger::GetInstance();
    pt.add("run.nFrames", prmMessenger->GetNumberOfFrames());
    pt.add("run.nParticles", prmMessenger->GetNumberOfParticles());
    
    //pt.add("geant4.random", fRM->GetRandomNumberStore());

    write_ini(fname, pt);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MpxDetector::UpdateSettings()
{
    G4RunManager *fRM = G4RunManager::GetRunManager();
    DetectorConstructionBase *myDet = (DetectorConstructionBase *)(fRM->GetUserDetectorConstruction());
    tpxMode = myDet->GetTpxMode();
}
