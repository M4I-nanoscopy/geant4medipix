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
/// \file ExportMgr.hh
/// \brief Definition of the ExportMgr class

#ifndef ExportMgr_h
#define ExportMgr_h 1

#include <list>
#include "ExportBase.hh"


#ifdef WITH_HDF5
#include "ExportHDF.hh"
#endif

#include "G4Types.hh"
#include "G4AutoLock.hh"
#include "MpxDetector.hh"

// class G4GenericMessenger;

class ExportMgr
{
public:
    /**
     * Adds pixel energy to HitsCollection
     * 
     */
    static ExportMgr *GetInstance();
    /**
     * constructor of the export manager class
     */
    ExportMgr();
    /**
     * The Export manager destructor
     */
    ~ExportMgr();
    /**
    * interaction data with sensor
    * \param HitsCollection* detector hits collection
    * \param event the eventID
    */
    void AddData(DetectorHitsCollection *, MpxDigitCollection *DigitCollection, G4int);
    /**
    * interaction data with sensor
    * \param HitsCollection* detector hits collection
    * \param event the eventID
    */
    void WriteData();
    /**
    * interaction data with sensor
    * \param 
    * 
    */
    void SetHDFFilename(G4String);

    void CreateDataFile();

private:
    static ExportMgr *instance;

    G4int lastEvent;/**< the last event */

    G4int nbEvents;/**< number of events */

    ExportBase *hdfExport; /**< the HDFExport instance */

    G4String filename; /**< the HDF5 filename*/
};


#endif
