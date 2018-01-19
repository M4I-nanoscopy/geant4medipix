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

#include "DetectorConstructionBase.hh"

#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"

#include <G4GenericMessenger.hh>
#include <PrimaryGeneratorAction.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExportHDF::ExportHDF()
{
    HitsCollectionCopy = new DetectorHitsCollection();
    DigitCollectionCopy = new MpxDigitCollection();
    filename    = "g4medipix.h5";
}

namespace
{
    G4Mutex HDF5Mutex = G4MUTEX_INITIALIZER;
    hid_t H5file;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExportHDF::AddSingleEvents(DetectorHitsCollection *HitsCollection)
{
    for (G4int i = 0; i < (G4int) HitsCollection->GetSize(); i++) {
        // Hard copy of object
        if (i > 0 && (*HitsCollection)[i]->GetTime() < (*HitsCollection)[i-1]->GetTime()) {
            return;
        }
        DetectorHit *hitDetector = (*HitsCollection) [i];
        DetectorHit *hitCopy = new DetectorHit(*hitDetector);

        HitsCollectionCopy->insert(hitCopy);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExportHDF::AddSingleDigits(MpxDigitCollection *DigitCollection)
{
    for (G4int i = 0; i < (G4int) DigitCollection->GetSize(); i++) {
        // Hard copy of object
        Digit *digit = (*DigitCollection) [i];
        DigitCollectionCopy->insert(new Digit(*digit));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExportHDF::AddEnergyPerPixel(DetectorHitsCollection *HitsCollection)
{
    G4int nbMatch = 0;
    G4int nbElements = (G4int) HitsCollection->GetSize();
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
            }
        }
    }

    delete[] matchArray;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportHDF::Write(G4String dataSetName) {
    size_t LENGTH = HitsCollectionCopy->GetSize();

    //G4cout << "Writing trajectories output per event to HDF5. Number of hits: " << LENGTH << G4endl;

    if ( LENGTH == 0) {
        return;
    }

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
    int exists = H5Lexists(file, dataSetName.c_str(), H5P_DEFAULT);

    if ( exists == 0 ) {
        H5Gcreate1(file, dataSetName.c_str(), sizeof(file));
    }

    DetectorConstructionBase *det = (DetectorConstructionBase *)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    if (! det->storeTraj) {
        return;
    }

    s1_t *s1 = new s1_t[LENGTH];

    G4int ev = (*HitsCollectionCopy)[0]->GetEvent();
    G4int temp = 0;

    for (size_t i = 0; i < LENGTH; i++) {

        DetectorHit *sensorHit = (*HitsCollectionCopy)[i];

        if (sensorHit->GetEvent() != ev || i == LENGTH - 1) {
            if (i == LENGTH - 1) {
                s1[temp].x = sensorHit->GetPosition().x() / nm;
                s1[temp].y = sensorHit->GetPosition().y() / nm;
                s1[temp].z = sensorHit->GetPosition().z() / nm;
                s1[temp].energy = sensorHit->GetEdep() / keV;
            }
            G4String tableName = dataSetName + std::to_string(ev);
            hsize_t dim[] = {(long long unsigned int) temp, 4};
            space = H5Screate_simple(2, dim, NULL);
            dataset = H5Dcreate(file, tableName, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dataset, H5T_IEEE_F64LE , H5S_ALL, H5S_ALL, H5P_DEFAULT, s1);
            H5Dclose(dataset);
            H5Sclose(space);

            ev = sensorHit->GetEvent();
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
    CloseOutputFile();
    HitsCollectionCopy = new DetectorHitsCollection();
}

void ExportHDF::WritePixels() {
    size_t number_digits = DigitCollectionCopy->GetSize();

    if (number_digits == 0) {
        return;
    }

    DetectorConstructionBase *det = (DetectorConstructionBase *)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    size_t nb = (size_t) det->GetNbPixels();

    // Reserve space
    G4double *pixels = (G4double*) calloc(2 * nb * nb, sizeof(G4double));
    // Set all ToA values to NaN and all ToT values to 0
    memset(pixels, -1, 2 * nb * nb * sizeof(G4double));
    memset(pixels, 0, 1 * nb * nb * sizeof(G4double));

    // Handles
    hid_t file = GetOutputFile();
    hid_t space, dataset;
    hsize_t dims[3] = {2, (hsize_t) nb, (hsize_t) nb};

    // Group pixels
    int exists = H5Lexists(file, "/g4medipix", H5P_DEFAULT);

    if ( exists == 0 ) {
        H5Gcreate1(file, "/g4medipix", sizeof(file));
    }

    // Get first event
    G4int event = (*DigitCollectionCopy)[0]->GetEvent();

    // Iterate over all entries in the list
    for (size_t i = 0; i < DigitCollectionCopy->GetSize(); i++) {

        Digit *d = (*DigitCollectionCopy)[i];

        // Digit collection start at 1 (hence -1)
        size_t x = (d->GetColumn() - 1)*nb + (d->GetLine() - 1);
        size_t y = nb*nb + (d->GetColumn() - 1)*nb + (d->GetLine() - 1);

        if ( event != d->GetEvent() || i == DigitCollectionCopy->GetSize() - 1) {

            if ( i == DigitCollectionCopy->GetSize() - 1 ) {
                pixels[x] = d->GetToT();
                pixels[y] = d->GetToA();
            }

            G4cout << "Writing sparse pixels output per event " << event << " to HDF5. Number of digits: " << i << G4endl;

            // Create dataset
            G4String tableName = "/g4medipix/" + std::to_string(event);
            space = H5Screate_simple(3, dims, NULL);
            dataset = H5Dcreate(file, tableName, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            // Write dataset
            H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pixels);

            H5Sclose(space);
            H5Dclose(dataset);

            // Set all ToA values to NaN and all ToT values to 0
            memset(pixels, -1, 2 * nb * nb * sizeof(G4double));
            memset(pixels, 0, 1 * nb * nb * sizeof(G4double));

            event = d->GetEvent();
        }

        pixels[x] = d->GetToT();
        pixels[y] = d->GetToA();
    }

    // Clean up
    free(pixels);
    delete DigitCollectionCopy;
    DigitCollectionCopy = new MpxDigitCollection();
    CloseOutputFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportHDF::SetFilename(G4String name)
{
// #ifdef G4MULTITHREADED
    // std::ostringstream os;
    // os << G4Threading::G4GetThreadId();
    //name.append("_t");
    //name.append(os.str());
// #endif

    filename = name;
}

void ExportHDF::SetAttributes() {
    hid_t file = GetOutputFile();

    DetectorConstructionBase *det = (DetectorConstructionBase *)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    PrimaryGeneratorAction *pga = (PrimaryGeneratorAction * )
            G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

    G4double energy =  pga->GetParticleGun()->GetParticleEnergy() / keV;
    G4double height = det->GetSensorThickness() / nm;
    G4String mat = det->GetSensorMaterial()->GetName();
    // Beam energy
    hid_t dataspace_id = H5Screate(H5S_SCALAR);
    hid_t att_energy = H5Acreate2 (file, "beam_energy", H5T_NATIVE_DOUBLE, dataspace_id,
                                   H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_energy,H5T_NATIVE_DOUBLE,&energy);
    // Sensor height
    hid_t att_height = H5Acreate2 (file, "sensor_height", H5T_NATIVE_DOUBLE, dataspace_id,
                                   H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_height,H5T_NATIVE_DOUBLE,&height);
    // Sensor material
    hid_t strr = H5Tcopy (H5T_C_S1);
    H5Tset_size (strr, strlen(mat));
    hid_t att_mat = H5Acreate2(file, "sensor_material", strr, dataspace_id,
                               H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite (att_mat, strr, mat);
    // Source
    H5Tset_size (strr, strlen("g4medipix"));
    hid_t att_source = H5Acreate2(file, "source", strr, dataspace_id,
                               H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite (att_source, strr, "g4medipix");

    // Clean up
    H5Sclose(dataspace_id);
    H5Aclose(att_energy);
    H5Aclose(att_height);
    H5Aclose(att_source);
    H5Aclose(att_mat);
    CloseOutputFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportHDF::CreateOutputFile() {
    G4cout << "Creating HDF5 output file " << filename.c_str() << G4endl;

    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    H5Fclose(file);
}

hid_t ExportHDF::GetOutputFile() {
    G4AutoLock lock(&HDF5Mutex);

    H5file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    return H5file;
}

void ExportHDF::CloseOutputFile() {
    H5Fflush(H5file, H5F_SCOPE_GLOBAL);
    H5Fclose(H5file);

    G4AutoLock unlock(&HDF5Mutex);

}

#endif
