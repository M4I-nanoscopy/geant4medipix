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
/// \file ExportRaw.hh
/// \brief Definition of the ExportRaw class


#include "DetectorHit.hh"
#include "ExportBase.hh"


#include <iostream>
#include <string>


/** ExportRaw class
* 
* this class is used to export HitsCollections to HDF5 files
*/
class ExportRaw : public ExportBase
{
public:
    G4String filename;
    /**
     * The ExportRaw constructor
     */
    ExportRaw();

    void AddSingleEvents(DetectorHitsCollection *);
    void AddSingleDigits(MpxDigitCollection *DigitCollection);
    /**
     *  Write trajectories
     * \param dataSetName the name of the dataset in the HDF5 file
     * \param event the event ID
     */
    void Write();
    /**
     *  Set hdf5 file name
     * \param name sets the name
     */
    void SetFilename(G4String);

    void CreateOutputFile();

    void WritePixels();

    void SetAttributes();

private:
    /** The hits collection copy from SD */
    DetectorHitsCollection *HitsCollectionCopy;
    MpxDigitCollection *DigitCollectionCopy;

    void CloseOutputFile();

    int GetOutputFile();
};

