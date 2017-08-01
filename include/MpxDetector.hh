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
/// \file MpxDetector.hh
/// \brief Definition of the MpxDetector class

#ifndef MpxDetector_h
#define MpxDetector_h 1

#include "Digit.hh"

#include "G4Types.hh"
#include "G4AutoLock.hh"
#include <list>
#include <inttypes.h>

#include "G4SystemOfUnits.hh"
class MpxDetector
{
public:
    static MpxDetector *GetInstance();

    ~MpxDetector();

    //interaction data with sensor
    void AddPixelEvents(MpxDigitCollection *);
    void WriteSparse();
    void WriteFrame();
    void ReInitMatrix();
    void ResetPixelMatrix();
    /**
     * Gets the new settings from the Detector Construction such as Threshold, Timepix mode etc.
     */
    void UpdateSettings();
    /**
     * Sets threshold dispersion matrix for the detector with random fluctuations
     */
    void SetThresholdDispersion();
    /**
     * Write the current simulation config to an ascii file
     */
    void WriteSimulationSettings();
    /**
     * Get the timepix threshold
     * \return threshold in keV
     */
    inline   G4double GetTpxThreshold() {
        return TpxThreshold / keV;
    }
    /**
     * Get the medipix threshold 1
     * \return threshold in keV
     */
    inline   G4double GetMpxThreshold1() {
        return MpxThreshold1 / keV;
    }
    /**
     * Get the medipix threshold 2
     * \return threshold in keV
     */
    inline   G4double GetMpxThreshold2() {
        return MpxThreshold2 / keV;
    }
    /**
     * Get the medipix threshold 3
     * \return threshold in keV
     */
    inline   G4double GetMpxThreshold3() {
        return MpxThreshold3 / keV;
    }
    /**
     * Get the medipix threshold 4
     * \return threshold in keV
     */
    inline   G4double GetMpxThreshold4() {
        return MpxThreshold4 / keV;
    }
    /**
     * set the timepix threshold
     * \param th in keV
     */
    inline void SetTpxThreshold(G4double th) {
        TpxThreshold = th;
    }
    /**
     * set the medipix threshold 1
     * \param th in keV
     */
    inline void SetMpxThreshold1(G4double th) {
        MpxThreshold1 = th;
    }
    /**
     * set the medipix threshold 2
     * \param th in keV
     */
    inline void SetMpxThreshold2(G4double th) {
        MpxThreshold2 = th;
    }
    /**
     * set the medipix threshold 3
     * \param th in keV
     */
    inline void SetMpxThreshold3(G4double th) {
        MpxThreshold3 = th;
    }
    /**
     * set the medipix threshold 4
     * \param th in keV
     */
    inline void SetMpxThreshold4(G4double th) {
        MpxThreshold4 = th;
    }

private:
    static MpxDetector *instance;

    MpxDetector();

    G4int writeCounter;
    
    G4int nPixel;/**< number of pixels nxn */
    G4int detectorType;
    G4bool csmMode;

    /** MP3RX counter for 4 thresholds
     * 
     */
    G4int *MpxPixelMat1;
    G4int *MpxPixelMat2;
    G4int *MpxPixelMat3;
    G4int *MpxPixelMat4;

    G4double MpxThreshold1;
    G4double MpxThreshold2;
    G4double MpxThreshold3;
    G4double MpxThreshold4;

    //Timepix
    G4double  TpxThreshold;
    G4double *TpxPixelMat;
    G4double *TpxThlDisp;
    G4int tpxMode;


    G4int frameCounter;

    /**
     * struct definition of the xyc+ output format
     * \param event
     * \param col x-coordinate of pixel
     * \param col y-coordinate of pixel
     * \param energy energy depotion(geant4)
     * \param tot time over threshold value after amplifier
     * \param toa time of arrival value from amplifier
     */
    struct snglEvent {
        uint32_t event;
        uint32_t col;
        uint32_t line;
        G4double energy;
        G4double tot;
        G4double toa;
    };
    
    std::list<snglEvent> sparseList;

};
#endif
