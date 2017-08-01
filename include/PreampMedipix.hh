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
/// \file include/PreampMedipix.hh
/// \brief Definition of the PreampMedipix class

#ifndef PreampMedipix_h
#define PreampMedipix_h 1

#include "MpxDetector.hh"
#include "boost/property_tree/ptree.hpp"
#include "DetectorConstructionBase.hh"

class PreampMedipix
{
public:
    /**
     * \brief the constructor
     * \param the boost ini-file tree
     */
    PreampMedipix(boost::property_tree::ptree pt);
    /**
     * \brief compute pixel response of induced charge signal
     * \param inducedPixelContent
     * \param event number
     * \return digitcollection of pixel events
     */
    MpxDigitCollection* GetPixelResponse(std::map<std::pair<G4int, G4int>, G4double * > *inducedPixelContent, G4int event);

    //*digitCollection
private:

    /**
     * pointer to MpxDetector
     */
    MpxDetector* detector;

    /**
     * Krummenacher current for leakage compensation (influences pulse shape)
     */
    G4double Ikrum;
    /**
     * stepsize of computation
     */
    G4double    pulsePrecision;
    /**
     * shutter time for induced pulse
     */
    G4double    maxPulseTime;
    /**
     * number of elements in induced charge array
     */
    G4int       nPulseArrayElements;
    /**
     * shutter for the preamp response
     */
    G4double    ampResponseTime;
    /**
     * number of elements in preamp output array
     */
    G4int       nAmpResponseElements;
    /**
     * number of elements in impulse response kernel
     */
    G4int       nImpResEl;
    /**
     * the impulse response kernel
     */
    G4double*   impRes;

    /**
     * number electron/hole pairs per energy
     */
    G4double    nElectronHolePairs;
    /**
     * lowest energy in output (cut very low energies)
     */
    G4double    minEnergy;

    /**
     * activate debug output
     */
    G4bool      DEBUG;
    /**
     * write induced and preamp pulses to files
     */
    G4bool      writePeakToFile;

    /**
     * the detector construction
     */
    DetectorConstructionBase* myDet;
    /**
     * activate charge summing mode
     */
    G4bool      useCSM;
    /**
     * choose preamp type. 0: charge integration, 1: convolution preamp
     */
    G4int       preampType;
    /**
     * amount of electronics noise in elementary charges
     */
    G4double elecSigma;

    /**
     * timepix threshold in keV used for ToT
     */
    G4double thresholdkeV;
    /**
     * timepix threshold in elementary charges used for ToT
     */
    G4double thresholdCharge;

    /**
     * threshold dispersion matrix
     */
    G4double* TpxThlDisp;

    /**
     * number of Pixels
     */
    G4int nPixel;
    
    /**
     * front end feedback capacitance
     */
    G4double Cf;

//    /**
//     * \brief timeOfArrival
//     */
//    G4double timeOfArrival;

//    /**
//     * \brief timeOverThreshold
//     */
//    G4double timeOverThreshold;

//    /**
//     * the maximal charge of the preamp output pulse
//     */
//    G4double peakCharge;

    /**
     * compute electronics noise for given charge
     * \param charge amount of charge
     */
    G4double ElectronicsNoise(G4double charge);

    /**
     * compute the convolution of input signal with preamp function. Returns the maximal Charge
     * \param in array of input signal
     * \param out array of the convolution output
     * \param length length of the input array
     * \param kernel the array of the kernel
     * \param kernel_length the length of the kernel
     */
    G4double convolve(G4double* in, G4double* out, G4int length, G4double* kernel, G4int kernel_length);
    /**
     * \brief compute ToT and peak charge of induced current signal
     * \param preamp
     * \param threshold
     * \return Tuple of ToT, PeakCharge and ToA
     */
    boost::tuple<G4double,G4double,G4double> getToTandPeak(G4double* preamp, G4double threshold, G4double thlDisp);

    /**
     * \brief PreampMedipix::SetTransferFunctions
     * sets the correct Transfer functions
     */
    void SetTransferFunctions();

    /**
     * \brief SetThresholdDispersion
     */
    void SetThresholdDispersion();


    /**
     * \brief integrate the induced charge and compute the energy
     * \param inducedPixelConten
     * \param event
     * \param digitCollection
     * \return the digitcollection is overwritten and can be returned to MpxDetector
     */
    void chargeIntegrationPreamp(std::map<std::pair<G4int, G4int>, G4double * > *inducedPixelConten, G4int event, MpxDigitCollection *digitCollection);
    /**
     * \brief convolute the induced charge with preamp response and compute the energy
     * \param inducedPixelConten
     * \param event
     * \param digitCollection
     * \return the digitcollection is overwritten and can be returned to MpxDetector
     */
    void convolutionPreamp(std::map<std::pair<G4int, G4int>, G4double * > *inducedPixelConten, G4int event, MpxDigitCollection *digitCollection);

    /**
     * \brief Mpx4like Preamp implementation
     * \param inducedPixelConten
     * \param event
     * \param digitCollection
     * \return the digitcollection is overwritten and can be returned to MpxDetector
     */
    void mpx4Preamp(std::map<std::pair<G4int, G4int>, G4double * > *inducedPixelConten, G4int event, MpxDigitCollection *digitCollection);

    /**
     * \brief write the peaks to file
     * \param filename
     * \param data array
     * \param nElements of data array
     */
    void writeToFile(G4String filename, G4double* data, G4int nElements);

};


#endif
