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
/// \file ExportHDF.hh
/// \brief Definition of the ExportHDF class

#ifndef ExportHDF_h
#define ExportHDF_h 1

#ifdef WITH_HDF5
#include "DetectorHit.hh"
#include "ExportBase.hh"
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
#endif

/** ExportHDF class
* 
* this class is used to export HitsCollections to HDF5 files
*/

class ExportHDF : public ExportBase
{
G4String filename;

public:
    /**
     * The ExportHDF constructor
     */
    ExportHDF();
    ExportHDF(G4String);

    /**
     * The ExportHDF destructor
     */
    virtual ~ExportHDF();
    /**
     * Export file function
     */
    void ExportToFile();
    
    /**
     * Adds pixel energy to HitsCollection
     * \param *HitsCollection DetectorHitsCollection to write to
     * \param event the eventID
     */
    void AddSingleEvents(DetectorHitsCollection *, G4int);
    /**
     *  Adds energy of single pixel to dataset
     * \param *HitsCollection DetectorHitsCollection to write to
     * \param event event to write
     */
    void AddEnergyPerPixel(DetectorHitsCollection *, G4int);
    /**
     *  Adds energy of single pixel to dataset
     * \param dataSetName the name of the dataset in the HDF5 file
     * \param event the event ID
     */
    void Write(G4String, G4int, G4double);
    /**
     * TODO
     */ 
    void WriteLast();
    /**
     *  Set hdf5 file name
     * \param name sets the name
     */
    void SetFilename(G4String);

private:
//     void DefineCommands();
    /** The hits collection copy from SD */
    DetectorHitsCollection *HitsCollectionCopy;
    /** */
    G4String    entryName;
    G4int       lastEvent;

    G4int       writeModulo;
    /** counter*/
    G4int       counter;

};

#endif
#endif
