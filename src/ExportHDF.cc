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

#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"

#include <typeinfo>
#include <sstream>
#include <G4GenericMessenger.hh>
#include <G4GeneralParticleSourceData.hh>
#include <G4GeneralParticleSourceMessenger.hh>
#include <G4GeneralParticleSource.hh>
#include <PrimaryGeneratorAction.hh>

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

void ExportHDF::Write(G4String dataSetName, G4int event) {
    size_t LENGTH = HitsCollectionCopy->GetSize();
    G4int RANK = 1;

    // structure  and dataset
    struct s1_t {
        G4int event;        //Event-Number
        G4int column;        //pixel
        G4int line;
        G4double x;        //position
        G4double y;
        G4double z;
        G4int particle;    //particle
        G4double tracklength;    //tracklength
        G4double energy;        //energy
    };

    // File datatype identifier
    hid_t s1_tid;
    // Handles
    hid_t file, dataset, space;

    // Get file
    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Group trajectories
    H5Gcreate1(file, "/trajectories", sizeof(file));

    //Attributes
    H5File h5File(filename.c_str(), H5F_ACC_RDWR);
    StrType str_type(PredType::C_S1, H5T_VARIABLE);
    DataSpace dspace(H5S_SCALAR);

    DetectorConstructionBase *det = (DetectorConstructionBase *)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    PrimaryGeneratorAction *pga = (PrimaryGeneratorAction * )
            G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

    if (!h5File.attrExists("beam_energy")) {
        Attribute att_energy = h5File.createAttribute("beam_energy", PredType::NATIVE_DOUBLE, dspace);
        G4double energy =  pga->GetParticleGun()->GetParticleEnergy() / keV;
        att_energy.write(PredType::NATIVE_DOUBLE, &energy);
    }
    if (!h5File.attrExists("sensor_height")) {
        G4double height = det->GetSensorThickness() * 1000000;
        Attribute att_height = h5File.createAttribute("sensor_height", PredType::NATIVE_DOUBLE, dspace);
        att_height.write(PredType::NATIVE_DOUBLE, &height);
    }
    if (!h5File.attrExists("sensor_material")) {
        G4String mat = det->GetSensorMaterial()->GetName();
        Attribute att_mat = h5File.createAttribute("sensor_material", str_type, dspace);
        att_mat.write(str_type, &mat);
    }

    s1_t *s1 = new s1_t[LENGTH];

    G4int ev = 0;
    G4int temp = 0;

    for (size_t i = 0; i < LENGTH; i++) {

        DetectorHit *sensorHit = (*HitsCollectionCopy)[i];

        if (sensorHit->GetEvent() != ev || i == LENGTH - 1) {
            if (i == LENGTH - 1) {
                s1[temp].event = sensorHit->GetEvent();
                s1[temp].column = sensorHit->GetColumn();
                s1[temp].line = sensorHit->GetLine();
                s1[temp].x = sensorHit->GetPosition().x();
                s1[temp].y = sensorHit->GetPosition().y();
                s1[temp].z = sensorHit->GetPosition().z();
                s1[temp].particle = sensorHit->GetParticleID();
                s1[temp].tracklength = sensorHit->GetTrackLength();
                s1[temp].energy = sensorHit->GetEdep() / keV;
                temp++;
            }
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
            G4String st = std::to_string(ev);
            G4String tableName = dataSetName + st;
            hsize_t dim[] = {(long long unsigned int) temp};
            space = H5Screate_simple(RANK, dim, NULL);
            dataset = H5Dcreate(file, tableName, s1_tid, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dataset, s1_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, s1);
            H5Dclose(dataset);
            H5Sclose(space);
            H5Tclose(s1_tid);

            ev++;
            temp = 0;
        }

        s1[temp].event = sensorHit->GetEvent();
        s1[temp].column = sensorHit->GetColumn();
        s1[temp].line = sensorHit->GetLine();
        s1[temp].x = sensorHit->GetPosition().x();
        s1[temp].y = sensorHit->GetPosition().y();
        s1[temp].z = sensorHit->GetPosition().z();
        s1[temp].particle = sensorHit->GetParticleID();
        s1[temp].tracklength = sensorHit->GetTrackLength();
        s1[temp].energy = sensorHit->GetEdep() / keV;
        temp++;
    }

    H5Fclose(file);

    delete HitsCollectionCopy;
    HitsCollectionCopy = new DetectorHitsCollection();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportHDF::SetFilename(G4String name)
{
    filename = name;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportHDF::WriteLast()
{
    Write(entryName, lastEvent);
}

#endif
