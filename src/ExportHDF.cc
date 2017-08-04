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
// $Id: ExportHDF.cc 71060 2013-06-10 15:03:19Z gcosmo $
//
/// \file src/ExportHDF.cc
/// \brief Implementation of the ExportHDF class

#ifdef WITH_HDF5

#include "ExportHDF.hh"

#include "RunAction.hh"
#include "DetectorConstructionBase.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"

#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif
#include <string>


#include <typeinfo>
#include <sstream>
#include <G4GenericMessenger.hh>
#include <G4GeneralParticleSourceData.hh>
#include <G4GeneralParticleSourceMessenger.hh>
#include <G4GeneralParticleSource.hh>
#include <G4ParticleGun.hh>

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExportHDF::ExportHDF()
{
    HitsCollectionCopy = new DetectorHitsCollection();
    lastEvent	= 0;
    //writeModulo = 1000;	//TODO now in Run.cc!!!
//     filename 	= "Medipix.h5";
    filename    = "Medipix.h5";
    entryName   = "trajectories";
    counter 	= 1;
    
//     DefineCommands();
}

ExportHDF::ExportHDF(G4String name)
{
    filename = name;
    HitsCollectionCopy = new DetectorHitsCollection();
    lastEvent   = 0;
    
    entryName   = "trajectories";
    counter     = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExportHDF::~ExportHDF() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExportHDF::AddSingleEvents(DetectorHitsCollection *HitsCollection, G4int event)
{
    for (G4int i = 0; i < (G4int) HitsCollection->GetSize(); i++) {

        //hard copy of object
        DetectorHit *hitDetector = (*HitsCollection) [i];
        DetectorHit *hitCopy = new DetectorHit(*hitDetector);

        HitsCollectionCopy->insert(hitCopy);
        counter++;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExportHDF::AddEnergyPerPixel(DetectorHitsCollection *HitsCollection, G4int event)
{
    G4int nbMatch = 0;
    G4int nbElements = HitsCollection->GetSize();
    G4int *matchArray = new G4int[nbElements];

    for (G4int i = 0; i < nbElements; i++) {
        matchArray[i] = 0;
    }

    while (nbMatch < nbElements) {

        for (G4int i = 0; i < nbElements; i++) {
            if (matchArray[i] == 0) {

                matchArray[i] = 1;
                nbMatch++;

                //new empty hit and fill with basic info
                DetectorHit *hitDetector = (*HitsCollection) [i];
                DetectorHit *hitCopyTmp = new DetectorHit();
                G4ThreeVector pos;
                hitCopyTmp->Add(hitDetector->GetEdep(), 0., pos, hitDetector->GetColumn(), hitDetector->GetLine(), hitDetector->GetEvent(),hitDetector->GetParticleID(),hitDetector->GetTime());

                for (G4int ii = i + 1; ii < nbElements; ii++) {
                    DetectorHit *hitDetectorCompare = (*HitsCollection) [ii];

                    if (hitDetectorCompare->GetLine() == hitCopyTmp->GetLine() && hitDetectorCompare->GetColumn() == hitCopyTmp->GetColumn()) {
                        //add energy in pixel
                        hitCopyTmp->SetEdep(hitCopyTmp->GetEdep() + hitDetectorCompare->GetEdep());
                        nbMatch++;
                        matchArray[ii] = 1;
                    }
                }
                HitsCollectionCopy->insert(hitCopyTmp);
                counter++;
            }
        }
    }

    delete[] matchArray;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportHDF::Write(G4String dataSetName, G4int event, G4double energy)
{
    G4int LENGTH = HitsCollectionCopy->GetSize();
    G4int RANK = 1;

    Exception::dontPrint();

    // structure  and dataset
    struct s1_t {
        G4int       event;		//Event-Number
        G4int       column;		//pixel
        G4int       line;
        G4double    x;		//position
        G4double    y;
        G4double    z;
        G4int       particle;	//particle
        G4double    tracklength;	//tracklength
        G4double    energy;		//energy
    };		//

    //typedef struct s1_t  s1_t;

    G4int sizes[event + 1];
    G4int ev;
    G4int comp=0;
    G4int cnter=0;
    //fill with data
    for (G4int i = 0; i < LENGTH; i++) {
        DetectorHit *sensorHit = (*HitsCollectionCopy)[i];
        ev = sensorHit->GetEvent();
        if (ev > comp) {
            sizes[comp] = cnter;
            comp = ev;
            cnter = 0;
        }
        cnter++;
    }
    sizes[event] = cnter;

    G4int sum_sizes[event + 1];
    for (G4int i = 0; i < event + 1; i++) {
        if (i==0) sum_sizes[i] = sizes[i];
        else sum_sizes[i] = sum_sizes[i-1] + sizes[i];
    }

    s1_t **s1;
    s1 = (s1_t**) malloc((event + 1)*sizeof(s1_t*));
    for (G4int i = 0; i < event + 1; i++) {
        s1[i] = (s1_t*) malloc(sizes[i]*sizeof(s1_t));
    }


    hid_t      s1_tid;     // File datatype identifier


    hid_t      file, dataset, space, faplist_id; /* Handles */


    //fill with data
    for (G4int i = 0; i < event + 1; i++) {
        for (G4int j = 0; j < sizes[i]; j++) {
            DetectorHit *sensorHit;
            if (i==0) sensorHit = (*HitsCollectionCopy)[j];
            else sensorHit = (*HitsCollectionCopy)[sum_sizes[i-1] + j];s1[i][j].event = sensorHit->GetEvent();
            s1[i][j].column = sensorHit->GetColumn();
            s1[i][j].line = sensorHit->GetLine();
            s1[i][j].x = sensorHit->GetPosition().x();
            s1[i][j].y = sensorHit->GetPosition().y();
            s1[i][j].z = sensorHit->GetPosition().z();
            s1[i][j].particle = sensorHit->GetParticleID();                            //FIXME
            s1[i][j].tracklength = sensorHit->GetTrackLength();                 //GetTime() / ns;
            s1[i][j].energy = sensorHit->GetEdep() / keV; //save in keV
        }
    }


    //get file
    try {
        faplist_id 	= H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_stdio(faplist_id);
        file 	= H5Fopen(filename.c_str(), H5F_ACC_RDWR, faplist_id);
    } catch (FileIException error) {
    }
    if (file < 0) {
        file 	= H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    H5Gcreate1(file,"/trajectories",sizeof(file));
    H5File file_( filename.c_str(), H5F_ACC_RDWR );

    StrType str_type(PredType::C_S1, H5T_VARIABLE);
    DataSpace dspace(H5S_SCALAR);
    DetectorConstructionBase* det = (DetectorConstructionBase*)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    if (!file_.attrExists("beam_energy") && energy != 0.0) {
        Attribute att_energy = file_.createAttribute("beam_energy",PredType::NATIVE_DOUBLE,dspace);
        G4double energy1 = energy*1000;
        att_energy.write(PredType::NATIVE_DOUBLE,&energy1);
    }
    if (!file_.attrExists("sensor_height")) {
        G4double height = det->GetSensorThickness()*1000000;
        Attribute att_height = file_.createAttribute("sensor_height",PredType::NATIVE_DOUBLE,dspace);
        att_height.write(PredType::NATIVE_DOUBLE,&height);
    }
    if (!file_.attrExists("sensor_material")) {
        G4String mat = det->GetSensorMaterial()->GetName();
        Attribute att_mat = file_.createAttribute("sensor_material",str_type,dspace);
        att_mat.write(str_type,&mat);
    }


    //create Datatype
    for (G4int i = 0; i < event + 1; i++) {
        //label conversion
        std::stringstream out;
        out << i;
        G4String currentEvent = out.str();
        G4String tableName = dataSetName + currentEvent;
        //herr_t     status;
        hsize_t    dim[] = { (long long unsigned int) sizes[i]};   /* Dataspace dimensions */
        //dataSpace
        space = H5Screate_simple(RANK, dim, NULL);

        s1_tid = H5Tcreate(H5T_COMPOUND, sizeof(s1_t));
        H5Tinsert(s1_tid, "Event", HOFFSET(s1_t, event), H5T_NATIVE_INT);
        H5Tinsert(s1_tid, "Column", HOFFSET(s1_t, column), H5T_NATIVE_INT);
        H5Tinsert(s1_tid, "Line", HOFFSET(s1_t, line), H5T_NATIVE_INT);
        H5Tinsert(s1_tid, "X", HOFFSET(s1_t, x), H5T_NATIVE_DOUBLE);
        H5Tinsert(s1_tid, "Y", HOFFSET(s1_t, y), H5T_NATIVE_DOUBLE);
        H5Tinsert(s1_tid, "Z", HOFFSET(s1_t, z), H5T_NATIVE_DOUBLE);
        H5Tinsert(s1_tid, "particle", HOFFSET(s1_t, particle), H5T_NATIVE_INT);
        H5Tinsert(s1_tid, "tracklength", HOFFSET(s1_t, tracklength), H5T_NATIVE_DOUBLE);
        H5Tinsert(s1_tid, "energy", HOFFSET(s1_t, energy), H5T_NATIVE_DOUBLE);

        //throws error when dataset already exists! maybe errorhandling.
        dataset = H5Dcreate(file, tableName, s1_tid, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataset, s1_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, s1[i]);

        H5Sclose(space);
        H5Dclose(dataset);

    }

    //close
    H5Tclose(s1_tid);
    H5Fclose(file);

    //empty list
    delete HitsCollectionCopy;
    HitsCollectionCopy 	= new DetectorHitsCollection();

    free(s1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportHDF::SetFilename(G4String name)
{
    filename = name;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportHDF::WriteLast()
{
    Write(entryName, lastEvent, 0.0);
}

#endif
