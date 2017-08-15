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

#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif

#include <string>

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

const int PIXELS_CHUNK_SIZE = 11;


/** ExportHDF class
* 
* this class is used to export HitsCollections to HDF5 files
*/
struct snglEvent {
    uint32_t event;
    uint32_t col;
    uint32_t line;
    G4double energy;
    G4double tot;
    G4double toa;
};
class ExportHDF : public ExportBase
{
public:
    G4String filename;
    /**
     * The ExportHDF constructor
     */
    ExportHDF();

    void AddSingleEvents(DetectorHitsCollection *);
    void AddSingleDigits(MpxDigitCollection *DigitCollection);
    /**
     *  Adds energy of single pixel to dataset
     * \param *HitsCollection DetectorHitsCollection to write to
     * \param event event to write
     */
    void AddEnergyPerPixel(DetectorHitsCollection *, G4int);
    /**
     *  Write trajectories
     * \param dataSetName the name of the dataset in the HDF5 file
     * \param event the event ID
     */
    void Write(G4String, G4int);
    /**
     *  Set hdf5 file name
     * \param name sets the name
     */
    void SetFilename(G4String);

    void CreateOutputFile();

    void WritePixels();

private:
//     void DefineCommands();
    /** The hits collection copy from SD */
    DetectorHitsCollection *HitsCollectionCopy;
    MpxDigitCollection *DigitCollectionCopy;
    /** */
    G4String    entryName;
    G4int       lastEvent;

    G4int       writeModulo;
    /** counter*/
    G4int       counter;

    G4int offset;

    void CloseOutputFile(hid_t file);

    hid_t GetOutputFile();

    std::list<snglEvent> sparseList;
};

#endif
#endif
