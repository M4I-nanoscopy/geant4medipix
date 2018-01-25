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
    filename    = "g4medipix.h5";
    MAX_TRAJ = 512;
}

namespace
{
    G4Mutex HDF5Mutex = G4MUTEX_INITIALIZER;
    hid_t H5file = -1;
    hid_t fapl_id;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExportHDF::AddSingleEvents(DetectorHitsCollection *HitsCollection)
{
    ExportHDF::Write(HitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExportHDF::AddSingleDigits(MpxDigitCollection *DigitCollection)
{
    ExportHDF::WritePixels(DigitCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportHDF::Write(DetectorHitsCollection *hc) {
    size_t LENGTH = hc->GetSize();

    if ( LENGTH == 0) {
        return;
    }

    // structure  and dataset
    struct s1_t {
        G4double x;
        G4double y;
        G4double z;
        G4double energy;
        G4double time;
    };

    // Handles
    hid_t file, dataset, trajSpace;

    auto *det = (DetectorConstructionBase *) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    if (! det->storeTraj) {
        return;
    }

    // Mutex lock
    G4AutoLock autoLock(&HDF5Mutex);

    // Get file
    file = GetOutputFile();

    // Array of structs to store
    auto *s1 = new s1_t[MAX_TRAJ];
    memset(s1, 0, MAX_TRAJ * sizeof(struct s1_t));

    G4int ev = (*hc)[0]->GetEvent();
    G4int temp = 0;

    for (size_t i = 0; i < LENGTH; i++) {

        if ( i == MAX_TRAJ ) {
            G4cout << "Found a trajectory larger than max trajectory size. Throwing it away" << G4endl;
            break;
        }

        DetectorHit *sensorHit = (*hc)[i];

        if (sensorHit->GetEvent() != ev || i == LENGTH - 1) {
            if (i == LENGTH - 1) {
                s1[temp].x = sensorHit->GetPosition().x() / nm;
                s1[temp].y = sensorHit->GetPosition().y() / nm;
                s1[temp].z = sensorHit->GetPosition().z() / nm;
                s1[temp].energy = sensorHit->GetEdep() / keV;
                s1[temp].time = sensorHit->GetTime() / picosecond;
            }

            G4String tableName = "/trajectories/" + std::to_string(ev);

            hsize_t dim[] = {(long long unsigned int) temp, 5};
            trajSpace = H5Screate_simple(2, dim, nullptr);
            dataset = H5Dcreate(file, tableName, H5T_IEEE_F64LE, trajSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            // G4cout << "Writing trajectories event " << ev << " output to HDF5. Number of hits: " << LENGTH << G4endl;

            H5Dwrite(dataset, H5T_IEEE_F64LE , H5S_ALL, H5S_ALL, H5P_DEFAULT, s1);
            H5Dclose(dataset);
            H5Sclose(trajSpace);

            memset(s1, 0, MAX_TRAJ * sizeof(struct s1_t));
            temp = 0;

            ev = sensorHit->GetEvent();

        }

        s1[temp].x = sensorHit->GetPosition().x() / nm;
        s1[temp].y = sensorHit->GetPosition().y() / nm;
        s1[temp].z = sensorHit->GetPosition().z() / nm;
        s1[temp].energy = sensorHit->GetEdep() / keV;
        s1[temp].time = sensorHit->GetTime() / picosecond;
        temp++;
    }

    // Close output (first close output, then destroy memory)
    CloseOutputFile();
    delete[] s1;
}

void ExportHDF::WritePixels(MpxDigitCollection *dc) {
    size_t number_digits = dc->GetSize();

    if (number_digits == 0) {
        return;
    }

    auto *det = (DetectorConstructionBase *)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    auto nb = (size_t) det->GetNbPixels();

    // Reserve space
    size_t pixel_size = 2 * nb * nb;
    auto *pixels = new G4double[pixel_size]{0};

    // Mutex lock
    G4AutoLock autoLock(&HDF5Mutex);

    // Handles
    hid_t file = GetOutputFile();
    hid_t dataset;

    // Get first event
    G4int event = (*dc)[0]->GetEvent();

    // Iterate over all entries in the list
    for (size_t i = 0; i < dc->GetSize(); i++) {

        Digit *d = (*dc)[i];

        // Used to have problems with eventID 0 popping up later. Not seen it for a while
        // if (d->GetEvent() == 0) {
        //    G4cout << "Found a digit with event 0. Throwing it away, not trusting it." << G4endl;
        //    continue;
        // }

        // Digit collection start at 1 (hence -1)
        size_t x = (d->GetColumn() - 1)*nb + (d->GetLine() - 1);
        size_t y = nb*nb + (d->GetColumn() - 1)*nb + (d->GetLine() - 1);

        if ( x > pixel_size || y > pixel_size) {
            G4cout << "Found a digit in event " << d->GetEvent() << " larger than pixel matrix. Throwing it away" << G4endl;
            continue;
        }

        if ( event != d->GetEvent() || i == dc->GetSize() - 1) {

            if ( i == dc->GetSize() - 1 ) {
                pixels[x] = d->GetToT();
                pixels[y] = d->GetToA();
            }

            // G4cout << "Writing sparse pixels output per event " << event << " to HDF5. Digits: " << i << G4endl;

            // Open dataset
            G4String tableName = "/g4medipix/" + std::to_string(event);
            dataset = H5Dopen1(file, tableName);

            // Write dataset
            H5Dwrite (dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pixels);
            H5Dclose(dataset);

            // Set all ToA values to 0 and all ToT values to 0
            memset(pixels, 0, pixel_size * sizeof(G4double));

            event = d->GetEvent();
        }

        pixels[x] = d->GetToT();
        pixels[y] = d->GetToA();
    }

    // Close output (first close output, then destroy memory)
    CloseOutputFile();
    delete[] pixels;
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
    // Mutex lock
    G4AutoLock autoLock(&HDF5Mutex);

    G4cout << "Setting attributes HDF5 output file " << filename.c_str() << G4endl;

    // Handles
    hid_t file = GetOutputFile();

    auto *det = (DetectorConstructionBase *) G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    G4double height = det->GetSensorThickness() / nm;
    G4double nbp = det->GetNbPixels();
    G4String mat = det->GetSensorMaterial()->GetName();
    hid_t dataspace_id = H5Screate(H5S_SCALAR);

    // Beam energy
    auto *pga = (PrimaryGeneratorAction * )G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
    G4double energy =  pga->GetParticleGun()->GetParticleEnergy() / keV;

    hid_t att_energy = H5Acreate2 (file, "beam_energy", H5T_NATIVE_DOUBLE, dataspace_id,
                               H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_energy,H5T_NATIVE_DOUBLE,&energy);

    // Sensor height
    hid_t att_height = H5Acreate2 (file, "sensor_height", H5T_NATIVE_DOUBLE, dataspace_id,
                                   H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_height,H5T_NATIVE_DOUBLE,&height);

    // Number of pixels
    hid_t nb_pixels = H5Acreate2 (file, "n_pixels", H5T_NATIVE_DOUBLE, dataspace_id,
                                   H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(nb_pixels,H5T_NATIVE_DOUBLE,&nbp);

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
    H5Aclose(nb_pixels);
    CloseOutputFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExportHDF::CreateOutputFile() {
    G4AutoLock autoLock(&HDF5Mutex);

    G4cout << "Creating HDF5 output file " << filename.c_str() << G4endl;

    hid_t clusterSpace, file, clusterId, clusterGroup, trajGroup;

    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl_id, H5F_CLOSE_SEMI);

    // Create file and dataset
    H5file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    clusterGroup = H5Gcreate1(H5file, "/g4medipix", sizeof(file));
    trajGroup = H5Gcreate1(H5file, "/trajectories", sizeof(file));

    // Get detector
    auto det = (DetectorConstructionBase *)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    auto nb = (size_t) det->GetNbPixels();

    // Number of events
    G4int n_events = G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed();

    // Mem space
    hsize_t dimsCluster[3] = {2, (hsize_t) nb, (hsize_t) nb};
    clusterSpace = H5Screate_simple(3, dimsCluster, nullptr);

    // Create dataset for each event
    for( G4int event = 0 ; event < n_events; event++) {
        G4String clusterName = "/g4medipix/" + std::to_string(event);
        clusterId = H5Dcreate(H5file, clusterName, H5T_IEEE_F64LE, clusterSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(clusterId);
    }

    H5Sclose(clusterSpace);
    H5Pclose(fapl_id);
    H5Gclose(clusterGroup);
    H5Gclose(trajGroup);

    H5Fflush(H5file, H5F_SCOPE_GLOBAL);
    H5Fclose(H5file);
    // G4cout << "Closed  " << H5file << " " << e << G4endl;
}

hid_t ExportHDF::GetOutputFile() {
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl_id, H5F_CLOSE_SEMI);

    H5file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl_id);


    return H5file;
}

void ExportHDF::CloseOutputFile() {
    H5Fflush(H5file, H5F_SCOPE_GLOBAL);

    H5Pclose(fapl_id);
    H5Fclose(H5file);

    // G4cout << "Closed  " << H5file << " " << e << G4endl;
}

#endif
