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
    filename    = "Medipix.h5";
    entryName   = "trajectories";
    counter 	= 1;
    offset = 0;
    
//     DefineCommands();
}

ExportHDF::ExportHDF(G4String name)
{
    filename = name;
    HitsCollectionCopy = new DetectorHitsCollection();
    lastEvent   = 0;
    
    entryName   = "trajectories";
    counter     = 1;
    offset = 0;
}

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

    // structure  and dataset
    struct s1_t {
        G4double x;
        G4double y;
        G4double z;
        G4double energy;
    };

    // Handles
    hid_t file, dataset, space;

    // Get file
    file = GetOutputFile();

    // Group trajectories
    H5Gcreate1(file, dataSetName.c_str(), sizeof(file));

    DetectorConstructionBase *det = (DetectorConstructionBase *)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    if (! det->storeTraj) {
        return;
    }

    s1_t *s1 = new s1_t[LENGTH];

    G4int ev = 0;
    G4int temp = 0;

    for (size_t i = 0; i < LENGTH; i++) {

        DetectorHit *sensorHit = (*HitsCollectionCopy)[i];

        if (sensorHit->GetEvent() != ev || i == LENGTH - 1) {
            if (i == LENGTH - 1) {
                s1[temp].x = sensorHit->GetPosition().x() / nm;
                s1[temp].y = sensorHit->GetPosition().y() / nm;
                s1[temp].z = sensorHit->GetPosition().z() / nm;
                s1[temp].energy = sensorHit->GetEdep() / keV;
                temp++;
            }
            G4String tableName = dataSetName + std::to_string(ev);
            hsize_t dim[] = {(long long unsigned int) temp, 4};
            space = H5Screate_simple(2, dim, NULL);
            dataset = H5Dcreate(file, tableName, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dataset, H5T_IEEE_F64LE , H5S_ALL, H5S_ALL, H5P_DEFAULT, s1);
            H5Dclose(dataset);
            H5Sclose(space);

            ev++;
            temp = 0;
        }

        s1[temp].x = sensorHit->GetPosition().x() / nm;
        s1[temp].y = sensorHit->GetPosition().y() / nm;
        s1[temp].z = sensorHit->GetPosition().z() / nm;
        s1[temp].energy = sensorHit->GetEdep() / keV;
        temp++;
    }

    delete[] s1;
    delete HitsCollectionCopy;
    CloseOutputFile(file);
    HitsCollectionCopy = new DetectorHitsCollection();
}

void ExportHDF::WritePixels(std::list<MpxDetector::snglEvent> list) {
    G4int number_events = G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed();
    if (number_events%PIXELS_CHUNK_SIZE == 0 && offset == number_events) return;

    DetectorConstructionBase *det = (DetectorConstructionBase *)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    G4int nb = det->GetNbPixels();

    // TODO: Get pixel size from settings
    G4double *pixels = (G4double*) calloc(PIXELS_CHUNK_SIZE * 2 * nb * nb, sizeof(G4double));

    // Handles
    hid_t file = GetOutputFile();
    hid_t data = PixelsDataset(file, number_events);

    G4int evv = offset;
    for (std::list<MpxDetector::snglEvent>::const_iterator iterator = list.begin(); iterator != list.end(); ++iterator) {
        struct MpxDetector::snglEvent e = *iterator;
        pixels[(e.event - offset)*2*nb*nb + e.col*nb + e.line] = e.tot;
        pixels[(e.event - offset)*2*nb*nb + nb*nb + e.col*nb + e.line] = e.toa;

        evv = e.event;
    }
    offset = evv + 1;

    hsize_t offset_[4];
    offset_[0] = ((offset-1)/PIXELS_CHUNK_SIZE)*PIXELS_CHUNK_SIZE;
    offset_[1] = 0;
    offset_[2] = 0;
    offset_[3] = 0;
    hsize_t size[4];
    size[0] = offset;
    size[1] = 2;
    size[2] = nb;
    size[3] = nb;
    H5Dset_extent(data,size);
    hid_t filespace = H5Dget_space (data);
    hsize_t dimsext[4];
    dimsext[0] = offset - ((offset-1)/PIXELS_CHUNK_SIZE)*PIXELS_CHUNK_SIZE;
    dimsext[1] = 2;
    dimsext[2] = nb;
    dimsext[3] = nb;
    H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset_, NULL,
                         dimsext, NULL);
    hid_t memspace = H5Screate_simple (4, dimsext, NULL);
    H5Dwrite (data, H5T_NATIVE_DOUBLE, memspace, filespace,
              H5P_DEFAULT, pixels);

    // Clean up
    free(pixels);
    H5Dclose (data);
    CloseOutputFile(file);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportHDF::SetFilename(G4String name)
{
    filename = name;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportHDF::CreateOutputFile() {

    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Attributes
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
        G4double height = det->GetSensorThickness() / nm;
        Attribute att_height = h5File.createAttribute("sensor_height", PredType::NATIVE_DOUBLE, dspace);
        att_height.write(PredType::NATIVE_DOUBLE, &height);
    }
    if (!h5File.attrExists("sensor_material")) {
        G4String mat = det->GetSensorMaterial()->GetName();
        Attribute att_mat = h5File.createAttribute("sensor_material", str_type, dspace);
        att_mat.write(str_type, &mat);
    }

    h5File.close();
    H5Fclose(file);
}

hid_t ExportHDF::GetOutputFile() {
    hid_t file;
    file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    return file;
}

void ExportHDF::CloseOutputFile(hid_t file) {
    H5Fclose(file);
}

hid_t ExportHDF::PixelsDataset(hid_t file, G4int nevents) {
    hid_t space, prop, dataset;
    DetectorConstructionBase *det = (DetectorConstructionBase *)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    G4int nb = det->GetNbPixels();
    G4int firstSize;

    /*if (nevents >= PIXELS_CHUNK_SIZE) {
        firstSize = PIXELS_CHUNK_SIZE;
    } else {
        firstSize = nevents;
    }*/

    hsize_t dims[4] = {(hsize_t) PIXELS_CHUNK_SIZE, 2, (hsize_t) nb, (hsize_t) nb};
    hsize_t maxDims[4] = {H5S_UNLIMITED, 2, (hsize_t) nb, (hsize_t) nb};

    space = H5Screate_simple (4, dims, maxDims);
    prop = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(prop, 4, dims);
    dataset = H5Dcreate2(file, "/pixels", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, prop, H5P_DEFAULT);
    H5Pclose(prop);
    H5Sclose(space);
    return dataset;
}

#endif
