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
// $Id: PreampMedipix.cc $
//
/// \file src/PreampMedipix.cc
/// \brief Implementation of the PreampMedipix class

#include "PreampMedipix.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"

#include <stdio.h>
#include <iostream>
#include <fstream>

#include <boost/concept_check.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

using namespace std;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PreampMedipix::PreampMedipix(boost::property_tree::ptree pt)
{

    G4RunManager *fRM = G4RunManager::GetRunManager();
    myDet = (DetectorConstructionBase *)(fRM->GetUserDetectorConstruction());
    useCSM = myDet->GetCsmMode();

    G4int detectorType = myDet->GetDetectorType();
    // set gain modes for Medipix3RX and preamp 
    if(detectorType == 0){
        if(useCSM == true){
            elecSigma = pt.get<G4double>("chip_medipix3rx.csm_enc");
        }else {
            elecSigma = pt.get<G4double>("chip_medipix3rx.spm_enc");
        }
        Cf = myDet->GetDetectorGain();
    } else if(detectorType == 1) { //Timepix
        elecSigma = pt.get<G4double>("chip_timepix.enc");
        Cf = 8e-15; // 8fF same as Medipix2MXR
    } else if(detectorType == 2) { // Dosepix
        elecSigma = pt.get<G4double>("chip_dosepix.enc");
        Cf = 8e-15;
    } else if (detectorType == 3) { // Timepix3
        elecSigma = pt.get<G4double>("chip_timepix3.enc");
        Cf = 3e-15; // 3fF
    }


    // get material and load material properties from ini file
    G4String material =  myDet->GetSensorMaterial()->GetName();
    G4String sensor;
    if (material == "G4_Si")
        sensor = "sensor_silicon";
    else if (material == "G4_CADMIUM_TELLURIDE")
        sensor = "sensor_cdte";

    pulsePrecision =    pt.get<G4double>("computation.pulsePrecision") * ns;  //must be ns!
    maxPulseTime =      pt.get<G4double>("computation.maxPulseTime") * ns;

    nPulseArrayElements = (G4int) maxPulseTime / pulsePrecision;
    ampResponseTime =   pt.get<G4double>("computation.ampResponseTime") * ns;
    nAmpResponseElements = (G4int)(ampResponseTime / pulsePrecision);

    nElectronHolePairs = pt.get<G4double>(sensor + ".nElectronHolePairs");
    minEnergy = 0.*keV;

    DEBUG = pt.get<G4bool>("testing.debug");
    writePeakToFile = pt.get<G4bool>("testing.writePeakToFile");

    preampType = pt.get<G4int>("computation.preampType");
    
    Ikrum = pt.get<G4double>("chip.Ikrum") * 1e-9;

    detector = MpxDetector::GetInstance();
    thresholdkeV = detector->GetTpxThreshold();
    thresholdCharge = thresholdkeV * 1000 / nElectronHolePairs ;

    //set the transfer functions
    SetTransferFunctions();

    //THL dispersion
    nPixel = myDet->GetNbPixels();
    TpxThlDisp = new G4double[nPixel * nPixel]();

    SetThresholdDispersion();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PreampMedipix::SetThresholdDispersion()
{
    //load .ini file with configuration data
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("DetectorConfig.ini", pt);

    //sensor properties
    G4String material =  myDet->GetSensorMaterial()->GetName();
    G4String sensor;
    if (material == "G4_Si")
        sensor = "sensor_silicon";
    else if (material == "G4_CADMIUM_TELLURIDE")
    {
        sensor = "sensor_cdte";
    }

    G4int detectorType = myDet->GetDetectorType();

    G4double disp = 0;

    if (detectorType == 0){ //Medipix3RX
        if(useCSM == true){
            disp = pt.get<G4double>("chip_medipix3rx.csm_threshold_dispersion");
        } else {
            disp = pt.get<G4double>("chip_medipix3rx.spm_threshold_dispersion");
        }

    }else if (detectorType == 1){ //Timepix1
        disp = pt.get<G4double>("chip_timepix.threshold_dispersion");

    }else if(detectorType == 2){ //Dosepix
        disp = pt.get<G4double>("chip_dosepix.threshold_dispersion");

    }else if (detectorType == 3) { // Timepix3
        disp = pt.get<G4double>("chip_timepix3.threshold_dispersion");
    }

    for (G4int i = 0; i < nPixel * nPixel; i++) {
        //G4cout << "DEBUG: Setting Threshold Dispersion " << disp << " i: " << i << " nPixel: " << nPixel*nPixel << G4endl;
        TpxThlDisp[i] = CLHEP::RandGauss::shoot(0, disp); //sigma of thl dispersion
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PreampMedipix::SetTransferFunctions(){
    
    
    //Preamp transfer function
    G4double Cd = 25e-15;  //Detector capacitance

    G4double Gmfb   = Ikrum/2.*20;  //Transconductance in the feedback loop (approximation from Weak Inversion)
    G4double Rf     = 2./Gmfb;       //feedback resistor

    G4double gm0=30e-6;  //transconductance parameter
    //G4double R0=20e6;    //Output resistance amplifier
    G4double C0=25e-15;  //Output capacitance amplifier
    G4double Ct=Cd*Cf+Cd*C0+Cf*C0;

    G4double rcpr=Ct/(gm0*Cf);
    G4double rcpf=Cf*Rf;

    G4double w1 = 1./rcpr;
    G4double w2 = 1./rcpf;


    //create transfer function:
    nImpResEl = nAmpResponseElements;
    impRes = new G4double[nImpResEl];
    for (G4int i = 0; i < nImpResEl; i++) {
        impRes[i] = exp(-i* pulsePrecision * 1e-9 *w2 - exp(-i* pulsePrecision * 1e-9* w1));
        //G4cout << "Debug1: " << i << " value: " << impRes[i] << G4endl;
    }
    //Shaper transfer function
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PreampMedipix::GetPixelResponse(map<pair<G4int, G4int>, G4double *> *inducedPixelContent, MpxDigitCollection* digitCollection, G4int event)
{
    thresholdkeV = detector->GetTpxThreshold();
    thresholdCharge = thresholdkeV * 1000 / nElectronHolePairs ;

    if(preampType == 0){
        chargeIntegrationPreamp(inducedPixelContent, event, digitCollection);
    } else if(preampType == 1){
        convolutionPreamp(inducedPixelContent, event, digitCollection);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PreampMedipix::chargeIntegrationPreamp(map<pair<G4int, G4int>, G4double *> *inducedPixelContent, G4int event, MpxDigitCollection *digitCollection)
{
    map<pair<G4int, G4int>, G4double> pixelsContent;

    map<pair<G4int, G4int>, G4double * >::iterator iPCItr;
    iPCItr = inducedPixelContent->begin();
    for (; iPCItr != inducedPixelContent->end() ; iPCItr++) {

        G4double maxCharge = 0;
        G4double *tempArray = (*iPCItr).second;
        G4double *ampResponse = new G4double[nAmpResponseElements]();

        //preamp response
        if(preampType == 0){
            for (G4int i = 0; i < nPulseArrayElements; i++) {
                maxCharge += tempArray[i];
            }
        }

        //convert maximum in preamp response back to energy //FIXME here TOT
        G4double energy = ElectronicsNoise(maxCharge) * nElectronHolePairs * eV;

        //safe new pixel energy, not yet Digit(SPM and CSM)
        if (energy > minEnergy) {
            std::pair<G4int, G4int> newPixel;
            newPixel.first  = (*iPCItr).first.first;
            newPixel.second = (*iPCItr).first.second;
            pixelsContent[newPixel] = energy;
        }
        //write preamp response to file
        if (writePeakToFile == true) {

            G4int k = (*iPCItr).first.first * 1000 + (*iPCItr).first.second;
            char filename[64];
            std::sprintf(filename, "preampSingle%d", k);

            writeToFile(filename, ampResponse, nAmpResponseElements);
        }

    }

    //print the energy in every pixel to the console
    if (DEBUG) {
        std::map<std::pair<G4int, G4int>, G4double>::iterator testItr =  pixelsContent.begin();
        G4cout << pixelsContent.size() << " number of entries in pixel Content" << G4endl;

        for (; testItr != pixelsContent.end() ; testItr++) {
            G4cout << (*testItr).first.first << " " << (*testItr).first.second << " " << (*testItr).second / keV << G4endl;
        }
        G4cout << G4endl;
    }

    //create digits. one per pixel
    std::map<std::pair<G4int, G4int>, G4double>::iterator pCItr =  pixelsContent.begin();
    std::map<std::pair<G4int, G4int>, G4double>::iterator tempPCItr =  pixelsContent.begin();

    //charge summing mode
    if (useCSM == true) {

        map<pair<G4int, G4int>, G4double> csmSumContent;
        std::pair<G4int, G4int> tempPixelCSM;
        std::pair<G4int, G4int> iteratePixel;

        //calculate the 2x2 sums and store in csmSumContent, equivalent of the sum/pixel in the chip
        for (; pCItr != pixelsContent.end() ; pCItr++) {

            G4double sumEnergy = 0;

            //get tempPixel
            tempPixelCSM = pCItr->first;
            iteratePixel = tempPixelCSM;

            //sumEnergy of pixel
            if (pixelsContent[tempPixelCSM] > 0) sumEnergy = pixelsContent[tempPixelCSM];

            iteratePixel.first += 1;
            tempPCItr = pixelsContent.find(iteratePixel);
            if (tempPCItr != pixelsContent.end()) {
                if (tempPCItr->second > 0) sumEnergy += tempPCItr->second;
            }

            iteratePixel.second -= 1;
            tempPCItr = pixelsContent.find(iteratePixel);
            if (tempPCItr != pixelsContent.end()) {
                if (tempPCItr->second > 0) sumEnergy += tempPCItr->second;
            }

            iteratePixel.first -= 1;
            tempPCItr = pixelsContent.find(iteratePixel);
            if (tempPCItr != pixelsContent.end()) {
                if (tempPCItr->second > 0) sumEnergy += tempPCItr->second;
            }

            //store the sum with the pixel ID
            csmSumContent[tempPixelCSM] = sumEnergy;
        }

        //compare the sums
        std::map<std::pair<G4int, G4int>, G4double>::iterator csmItr =  csmSumContent.begin();
        for (; csmItr != csmSumContent.end() ; csmItr++) {

            //get tempPixel
            tempPixelCSM = csmItr->first;
            iteratePixel = tempPixelCSM;

            //neighbors
            G4double cmpEnergy1 = 0;
            G4double cmpEnergy2 = 0;
            G4double cmpEnergy3 = 0;
            G4double cmpEnergy4 = 0;

            //diagonal
            G4double cmpEnergy5 = 0;
            G4double cmpEnergy6 = 0;
            G4double cmpEnergy7 = 0;
            G4double cmpEnergy8 = 0;

            //get sumEnergy of pixel
            G4double pixelEnergy = csmSumContent[tempPixelCSM];

            iteratePixel = tempPixelCSM;
            iteratePixel.first += 1;
            tempPCItr = csmSumContent.find(iteratePixel);
            if (tempPCItr != csmSumContent.end()) {
                cmpEnergy1 += tempPCItr->second;
            }

            iteratePixel = tempPixelCSM;
            iteratePixel.first -= 1;
            tempPCItr = csmSumContent.find(iteratePixel);
            if (tempPCItr != csmSumContent.end()) {
                cmpEnergy2 += tempPCItr->second;
            }

            iteratePixel = tempPixelCSM;
            iteratePixel.second += 1;
            tempPCItr = csmSumContent.find(iteratePixel);
            if (tempPCItr != csmSumContent.end()) {
                cmpEnergy3 += tempPCItr->second;
            }

            iteratePixel = tempPixelCSM;
            iteratePixel.second -= 1;
            tempPCItr = csmSumContent.find(iteratePixel);
            if (tempPCItr != csmSumContent.end()) {
                cmpEnergy4 += tempPCItr->second;
            }

            iteratePixel = tempPixelCSM;
            iteratePixel.first -= 1;
            iteratePixel.second -= 1;
            tempPCItr = csmSumContent.find(iteratePixel);
            if (tempPCItr != csmSumContent.end()) {
                cmpEnergy5 += tempPCItr->second;
            }

            iteratePixel = tempPixelCSM;
            iteratePixel.first += 1;
            iteratePixel.second -= 1;
            tempPCItr = csmSumContent.find(iteratePixel);
            if (tempPCItr != csmSumContent.end()) {
                cmpEnergy6 += tempPCItr->second;
            }

            iteratePixel = tempPixelCSM;
            iteratePixel.first -= 1;
            iteratePixel.second += 1;
            tempPCItr = csmSumContent.find(iteratePixel);
            if (tempPCItr != csmSumContent.end()) {
                cmpEnergy7 += tempPCItr->second;
            }

            iteratePixel = tempPixelCSM;
            iteratePixel.first += 1;
            iteratePixel.second += 1;
            tempPCItr = csmSumContent.find(iteratePixel);
            if (tempPCItr != csmSumContent.end()) {
                cmpEnergy8 += tempPCItr->second;
            }

            //create digit
            if (pixelEnergy > cmpEnergy1 && pixelEnergy > cmpEnergy2 && pixelEnergy > cmpEnergy3 && pixelEnergy > cmpEnergy4 && pixelEnergy > cmpEnergy5 && pixelEnergy > cmpEnergy6 && pixelEnergy > cmpEnergy7 && pixelEnergy > cmpEnergy8) {
                Digit *digit = new Digit;
                digit->SetColumn(tempPixelCSM.first);
                digit->SetLine(tempPixelCSM.second);
                digit->SetEnergy(pixelEnergy / keV);
                digit->SetEvent(event);
                digitCollection->insert(digit);
            }
        }

    } else {  //no CSM
        for (; pCItr != pixelsContent.end() ; pCItr++) {
            //convert all events into digits
            G4double energy = (*pCItr).second;
            if (energy > minEnergy) {
                Digit *digit = new Digit;
                digit->SetColumn((*pCItr).first.first);
                digit->SetLine((*pCItr).first.second);
                digit->SetEnergy((energy / keV));
                digit->SetEvent(event);
                digitCollection->insert(digit);
            }
        }
    }

    //free memory
    iPCItr = inducedPixelContent->begin();
    for (; iPCItr != inducedPixelContent->end() ; iPCItr++) {
        G4double *tempArray = (*iPCItr).second;
        delete[] tempArray;
    }
    inducedPixelContent->clear();
    delete inducedPixelContent;
    pixelsContent.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PreampMedipix::convolutionPreamp(map<pair<G4int, G4int>, G4double *> *inducedPixelContent, G4int event, MpxDigitCollection *digitCollection)
{
    //SPM mode
    if(useCSM == false){

        map<pair<G4int, G4int>, G4double * >::iterator iPCItr;
        iPCItr = inducedPixelContent->begin();

        //iterate over induced currents in pixels
        for (; iPCItr != inducedPixelContent->end() ; iPCItr++) {

            G4double *ampResponse = new G4double[nAmpResponseElements]();
            G4double *inducedChargeArray = (*iPCItr).second;

            //Preamp
            convolve(inducedChargeArray, ampResponse, nPulseArrayElements, impRes, nImpResEl);

            //noise
            for(G4int i=0; i<nAmpResponseElements; i++){
                ampResponse[i] = ElectronicsNoise(ampResponse[i]);
            }

            //getEnergy
            boost::tuple<G4double,G4double,G4double> totPeakToa;
            G4int pixelIndex = ((*iPCItr).first.first) * nPixel + ((*iPCItr).first.second);
            G4double thlDisp = TpxThlDisp[pixelIndex];

            totPeakToa = getToTandPeak(ampResponse, thresholdCharge, thlDisp);
            delete[] ampResponse;

            //convert to digits from peak value //FIXME here ToT
            if (boost::get<1>(totPeakToa) > minEnergy) {
                Digit *digit = new Digit;
                digit->SetColumn((*iPCItr).first.first);
                digit->SetLine((*iPCItr).first.second);
                digit->SetEnergy((boost::get<1>(totPeakToa) / keV));
                digit->SetToT(boost::get<0>(totPeakToa));
                digit->SetToA(boost::get<2>(totPeakToa));
                digit->SetEvent(event);
                digitCollection->insert(digit);
            }
        }

        //free memory
        iPCItr = inducedPixelContent->begin();
        for (; iPCItr != inducedPixelContent->end() ; iPCItr++) {
            G4double *tempArray = (*iPCItr).second;
            delete[] tempArray;
        }
        inducedPixelContent->clear();
        delete inducedPixelContent;
    }

    //CSM mode
    if(useCSM == true){

        map<pair<G4int, G4int>, G4double * >::iterator iPCItr;
        iPCItr = inducedPixelContent->begin();

        map<pair<G4int, G4int>, G4double * >* preampResponse    = new map<pair<G4int, G4int>, G4double * >;
        map<pair<G4int, G4int>, boost::tuple<G4double,G4double,G4double> >* sumPreampResponse = new map<pair<G4int, G4int>, boost::tuple<G4double,G4double,G4double> >;
        pair<G4int,G4int> tempPixel;
        G4double* emptyDummy = new G4double[nAmpResponseElements]();

        //preampresponse of every pixel
        for (; iPCItr != inducedPixelContent->end() ; iPCItr++) {
            G4double* ampResponse = new G4double[nAmpResponseElements]();
            G4double* inducedChargeArray = (*iPCItr).second;

            //Preamp
            convolve(inducedChargeArray, ampResponse, nPulseArrayElements, impRes, nImpResEl);

            //write single preamp response to file
            if (writePeakToFile == true) {

                G4int k = (*iPCItr).first.first * 1000 + (*iPCItr).first.second;
                char filename[64];
                std::sprintf(filename, "preampSingle%d", k);

                writeToFile(filename, ampResponse, nAmpResponseElements);
            }

            //store single preampResponse
            tempPixel.first = (*iPCItr).first.first;
            tempPixel.second = (*iPCItr).first.second;
            preampResponse->operator [](tempPixel) = ampResponse;
        }

        //sum the responses and store
        map<pair<G4int, G4int>, G4double * >::iterator pRItr;
        map<pair<G4int, G4int>, G4double * >::iterator tempPRItr;
        pRItr = preampResponse->begin();

        for(; pRItr != preampResponse->end(); pRItr++){

            pair<G4int,G4int>mainPixel = pRItr->first;
            pair<G4int,G4int>iteratePixel = mainPixel;

            //get the preamp arrays
            G4double *tempPreampReponse;
            //sum the responses
            G4double* csmSumPreamp = new G4double[nAmpResponseElements]();

            G4int sumArrLow = 0;
            G4int sumArrHigh = 1;

            for(G4int i=sumArrLow; i<=sumArrHigh; i++){
                for(G4int j=sumArrLow; j<=sumArrHigh; j++){
                    //set pixel
                    iteratePixel.first = mainPixel.first + i;
                    iteratePixel.second = mainPixel.second + j;

                    //find response pulse
                    tempPRItr = preampResponse->find(iteratePixel);
                    if(tempPRItr != preampResponse->end()){
                        tempPreampReponse = tempPRItr->second;

                        //add to sum
                        for(G4int k=0; k < nAmpResponseElements; k++){
                            csmSumPreamp[k] += tempPreampReponse[k];
                        }
                    }
                }
            }

            //add noise
            for(G4int i=0; i < nAmpResponseElements; i++){
                csmSumPreamp[i] = ElectronicsNoise(csmSumPreamp[i]);
            }

            //write CSM preamp response to file
            if (writePeakToFile == true) {

                G4int k = (*pRItr).first.first * 1000 + (*pRItr).first.second;
                char filename[64];
                std::sprintf(filename, "preampCSM%d", k);

                writeToFile(filename, csmSumPreamp, nAmpResponseElements);
            }

            //get energy
            boost::tuple<G4double,G4double,G4double> totPeakToa;
            G4int pixelIndex = ((*pRItr).first.first) * nPixel + ((*pRItr).first.second);
            G4double thlDisp = TpxThlDisp[pixelIndex];

            totPeakToa = getToTandPeak(csmSumPreamp, thresholdCharge, thlDisp);

            sumPreampResponse->operator [](mainPixel) = totPeakToa;

            //delete the array
            delete[] csmSumPreamp;
        }

        //compare the sums
        map<pair<G4int, G4int>, boost::tuple<G4double,G4double,G4double> >::iterator sPRItr = sumPreampResponse->begin();
        map<pair<G4int, G4int>, boost::tuple<G4double,G4double,G4double> >::iterator tempSPRItr;
        for (; sPRItr != sumPreampResponse->end() ; sPRItr++) {

            //get tempPixel
            pair<G4int,G4int> tempPixelCSM = sPRItr->first;
            pair<G4int,G4int> iteratePixel = tempPixelCSM;

            //first value: ToT, second: peak amplitude
            boost::tuple<G4double,G4double,G4double> cmpEnergy(0,0,0);


            //get sumEnergy of pixel
            boost::tuple<G4double,G4double,G4double> pixelEnergy = sumPreampResponse->operator [](tempPixelCSM);
            G4bool isMax = true;

            //sets the borders for the comparison
            G4int compArrLow    = -1;
            G4int compArrHigh   = 1;


            //compare neighbouring pixels with main pixel
            for(G4int i=compArrLow; i<=compArrHigh; i++){
                for(G4int j=compArrLow; j<=compArrHigh; j++){
                    if(!(i == 0 && j == 0)){
                        iteratePixel = tempPixelCSM;
                        iteratePixel.first += i;
                        iteratePixel.second += j;
                        tempSPRItr = sumPreampResponse->find(iteratePixel);
                        if (tempSPRItr != sumPreampResponse->end()) {
                            cmpEnergy = tempSPRItr->second;
                            if(boost::get<1>(cmpEnergy) > boost::get<1>(pixelEnergy)) isMax = false;
                        }
                    }
                }
            }

            //create digit
            if (isMax == true) {
                Digit *digit = new Digit;
                digit->SetColumn(tempPixelCSM.first);
                digit->SetLine(tempPixelCSM.second);
                digit->SetEnergy((boost::get<1>(pixelEnergy) / keV));
                digit->SetToT(boost::get<0>(pixelEnergy));
                digit->SetToA(boost::get<2>(pixelEnergy));
                digit->SetEvent(event);
                digitCollection->insert(digit);
            }
        }

        //free memory
        //induced pulses
        iPCItr = inducedPixelContent->begin();
        for (; iPCItr != inducedPixelContent->end() ; iPCItr++) {
            G4double *tempArray = (*iPCItr).second;
            delete[] tempArray;
        }
        inducedPixelContent->clear();
        delete inducedPixelContent;


        //preamp response
        pRItr = preampResponse->begin();

        for(; pRItr != preampResponse->end(); pRItr++){
            G4double* tempArray = (*pRItr).second;
            delete[] tempArray;
        }
        preampResponse->clear();
        delete preampResponse;

        //preamp sum
        sumPreampResponse->clear();
        delete sumPreampResponse;

        delete[] emptyDummy;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double PreampMedipix::ElectronicsNoise(G4double charge)
{
    G4double sigma = elecSigma;
    if(useCSM == true && preampType != 2) sigma = sigma * 2;
    return CLHEP::RandGauss::shoot(charge, sigma);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double PreampMedipix::convolve(G4double* in, G4double* out, G4int inLength, G4double* kernel, G4int kernel_length)
{
    G4double maxCharge = 0;

    for(G4int i=0; i<(kernel_length); i++){

        out[i] = 0.0;
        for(G4int k=0; k<kernel_length; k++){
            if((i-k) >=0 && (i-k) < inLength){
                out[i] += in[i-k] * kernel[k];
            }
        }
        if(out[i] > maxCharge) maxCharge = out[i];
    }

    return maxCharge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
boost::tuple<G4double,G4double,G4double> PreampMedipix::getToTandPeak(G4double* preamp, G4double threshold, G4double thlDisp){
    G4int risingEdge = 0;

    boost::tuple<G4double,G4double,G4double> returnValues(0,0,0);
    // 0 = ToT, 1 = peakCharge, 2 = ToA

    for(G4int i=0; i<nAmpResponseElements; i++){

        if(risingEdge == 0 && preamp[i] >= (threshold+thlDisp)){
            risingEdge = i;
            // Set timeOfArrival
            if(boost::get<2>(returnValues) == 0) {
                boost::get<2>(returnValues) = i * pulsePrecision;
            }
        }

        if(risingEdge != 0 && preamp[i] <= (threshold+thlDisp)){
            // Set timeOverThreshold
            if((i - risingEdge) * pulsePrecision > boost::get<0>(returnValues)) {
                boost::get<0>(returnValues) = (i - risingEdge) * pulsePrecision;
            }
            risingEdge = 0;
        }

        // set the peakCharge
        if(boost::get<1>(returnValues) < (preamp[i]+thlDisp)) boost::get<1>(returnValues) = (preamp[i]+thlDisp);
    }

    boost::get<1>(returnValues) = boost::get<1>(returnValues) * nElectronHolePairs * eV;

    return returnValues;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PreampMedipix::writeToFile(G4String filename, G4double* data, G4int nElements){
    ofstream myfile;

    myfile.open(filename);

    for (G4int it = 0; it < nElements; it++) {
        myfile << data[it] << " ";
    }
    myfile.close();
}
