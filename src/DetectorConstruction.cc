////
//// ********************************************************************
//// * License and Disclaimer                                           *
//// *                                                                  *
//// * The  Geant4 software  is  copyright of the Copyright Holders  of *
//// * the Geant4 Collaboration.  It is provided  under  the terms  and *
//// * conditions of the Geant4 Software License,  included in the file *
//// * LICENSE and available at  http://cern.ch/geant4/license .  These *
//// * include a list of copyright holders.                             *
//// *                                                                  *
//// * Neither the authors of this software system, nor their employing *
//// * institutes,nor the agencies providing financial support for this *
//// * work  make  any representation or  warranty, express or implied, *
//// * regarding  this  software system or assume any liability for its *
//// * use.  Please see the license in the file  LICENSE  and URL above *
//// * for the full disclaimer and the limitation of liability.         *
//// *                                                                  *
//// * This  code  implementation is the result of  the  scientific and *
//// * technical work of the GEANT4 collaboration.                      *
//// * By using,  copying,  modifying or  distributing the software (or *
//// * any work based  on the software)  you  agree  to acknowledge its *
//// * use  in  resulting  scientific  publications,  and indicate your *
//// * acceptance of all terms of the Geant4 Software license.          *
//// ********************************************************************
////
//// $Id: DetectorConstruction.cc 70424 2013-05-30 09:11:36Z gcosmo $
////
///// \file src/DetectorConstruction.cc
///// \brief Implementation of the DetectorConstruction class

//#include "DetectorConstruction.hh"
//#include "G4SystemOfUnits.hh"

//#include "G4Box.hh"
//#include "G4LogicalVolume.hh"
//#include "G4Tubs.hh"
//#include "G4SubtractionSolid.hh"
//#include "G4PVPlacement.hh"
//#include "G4PVReplica.hh"
//#include "G4AssemblyVolume.hh"
//#include "G4Colour.hh"

//// CADMESH //
//#define NOVCGLIB
//#include "CADMesh.hh"

//#include "G4NistManager.hh"
//#include "G4VisAttributes.hh"
//#include "G4PhysicalConstants.hh"
//#include "G4SystemOfUnits.hh"

////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//G4VPhysicalVolume *DetectorConstruction::DefineVolumes()
//{
//    CalculateGeometry();

//    G4double worldSizeXY = 1.2 * sensorSizeXY;
//    G4double worldSizeZ = 10 * cm;

//    // Get materials
//    G4Material *defaultMaterial = G4Material::GetMaterial("G4_Galactic");

//    if (! defaultMaterial || ! sensorMaterial) {
//        G4cerr << "Cannot retrieve materials already defined. " << G4endl;
//        G4cerr << "Exiting application " << G4endl;
//        exit(1);
//    }

//    //
//    // World
//    //
//    G4VSolid *worldS
//        = new G4Box("World",           // its name
//                    worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2); // its size

//    G4LogicalVolume *worldLV = new G4LogicalVolume(worldS,         // its solid
//            worldMaterial,// its material
//            "World");       // its name

//    G4VPhysicalVolume *worldPV = new G4PVPlacement(0,        // no rotation
//            G4ThreeVector(),  // at (0,0,0)
//            worldLV,          // its logical volume
//            "World",          // its name
//            0,                // its mother  volume
//            false,            // no boolean operation
//            0,                // copy number
//            fCheckOverlaps);  // checking overlaps

//// The Sensor:
//    //pixel
////    G4Box *pixelSolid = new G4Box("pixel_cell",
////                                  pixelSize / 2.,
////                                  pixelSize / 2.,
////                                  sensorThickness / 2);

//    G4ThreeVector offset = G4ThreeVector(0, 0, 0);
//    CADMesh * mesh = new CADMesh("test.stl", "STL", mm, offset, false);

//    G4VSolid* pixelSolid = mesh->TessellatedMesh();

//    G4LogicalVolume *pixelLogicalVolume = new G4LogicalVolume(pixelSolid,
//            sensorMaterial,
//            "pixel_cell");

//    // row
//    G4Box *rowSolid = new G4Box("row",
//                                (pixelSize * nPixel) / 2.,
//                                pixelSize / 2.,
//                                sensorThickness / 2.);

//    G4LogicalVolume *rowLogicalVolume = new G4LogicalVolume(rowSolid,
//            defaultMaterial,
//            "row");

//    new G4PVReplica("pixel_cell",
//                    pixelLogicalVolume,
//                    rowLogicalVolume,
//                    kXAxis,
//                    nPixel,
//                    pixelSize);

//    //sensor layer
//    G4Box *sensorSolid = new G4Box("Sensor",
//                                   (pixelSize * nPixel) / 2.,
//                                   (pixelSize * nPixel) / 2.,
//                                   sensorThickness / 2.);

//    G4LogicalVolume *sensorLogicalVolume = new G4LogicalVolume(sensorSolid,
//            defaultMaterial,
//            "Sensor");

//    new G4PVReplica("row",
//                    rowLogicalVolume,
//                    sensorLogicalVolume,
//                    kYAxis,
//                    nPixel,
//                    pixelSize);

//    new G4PVPlacement(G4Transform3D(fRotation,            //no rotation
//                                    G4ThreeVector(0, 0, 0)), //at (0,0,0)
//                      sensorLogicalVolume,      //its logical volume
//                      "sensor",      //its name
//                      worldLV,       //its mother  volume
//                      false,         //no boolean operation
//                      0);            //copy number



//    if (filter == true) {
//        //Fileters go here
//        G4Material *filterMaterial1 = G4Material::GetMaterial("G4_Al");
//        G4Box *filter1 = new G4Box("filterBox", (pixelSize * nPixel) / 2., (pixelSize * nPixel) / 2., filterThickness / 2.);
//        G4LogicalVolume *filterLV1 = new G4LogicalVolume(filter1,
//                filterMaterial,
//                "filLV");

//        // the translation translation vector
//        G4ThreeVector tFilter(0, 0, -(filterZ + filterThickness / 2 + sensorThickness / 2));

//        new G4PVPlacement(G4Transform3D(fRotation, tFilter),
//                          filterLV1,
//                          "filter1_",
//                          worldLV,
//                          false,
//                          0);
//    }

////-----------------------------------------------------------------Collimator
////--------------------------------------------------------------------------
//    if (collimator == true) {

//        G4Box  *colBox = new G4Box("collimator_box",
//                                   pixelSize / 2.,
//                                   pixelSize / 2.,
//                                   collimatorThickness / 2);

//        G4Tubs *colCyl =
//            new G4Tubs("collimator_cylinder", 0, colRadii, collimatorThickness, 0, twopi); // r:     0 mm -> 50 mm
//        // z:   -50 mm -> 50 mm
//        // phi:   0 ->  2 pi

//        G4SubtractionSolid *collimatorElementSolid =
//            new G4SubtractionSolid("collimator_cell", colBox, colCyl);


//        G4LogicalVolume *collimatorElementLogigalVolume = new G4LogicalVolume(collimatorElementSolid,
//                collimatorMaterial,
//                "collimator_cell");
//        G4AssemblyVolume *assemblyCollimator =  new G4AssemblyVolume();

//        G4RotationMatrix *Rcollimator = new G4RotationMatrix();
//        G4ThreeVector *TcollimatorE = new G4ThreeVector(0, 0, -1000 * um);
//        // make one assembly out of all collimators
//        for (int i = 0; i < nPixel; i++) {
//            for (int j = 0; j < nPixel; j++) {
//                TcollimatorE->setX(-pixelSize * nPixel / 2.0 + pixelSize / 2.0 + i * pixelSize);
//                TcollimatorE->setY(-pixelSize * nPixel / 2.0 + pixelSize / 2.0 + j * pixelSize);

//                assemblyCollimator->AddPlacedVolume(collimatorElementLogigalVolume, *TcollimatorE, Rcollimator);
//            }
//        }
//        G4ThreeVector *Tcollimator = new G4ThreeVector();
//        assemblyCollimator->MakeImprint(worldLV, *Tcollimator, &fRotation);

//    }//end if(collimator == true)

//    //-----------------------------------------------------------------Bumps
//    //--------------------------------------------------------------------------
//    if (bumps == true) {
//        G4Tubs *bumpCyl = new G4Tubs("bump_cyl", 0, bumpRadii, bumpHeight / 2.0, 0, twopi);

//        G4LogicalVolume *bumpLogicalVolume = new G4LogicalVolume(bumpCyl,
//                bumpMaterial,
//                "bump_cell");
//        // make one assembly out of all bumps
//        G4AssemblyVolume *assemblyBumps =  new G4AssemblyVolume();

//        G4RotationMatrix *Rbump = new G4RotationMatrix();
//        G4ThreeVector *Tbump = new G4ThreeVector();
//        Tbump->setZ(sensorThickness / 2.0 + bumpHeight / 2.0); // z stays constant
//        // add every bump to assembly
//        for (int i = 0; i < nPixel; i++) {
//            for (int j = 0; j < nPixel; j++) {
//                Tbump->setX(-pixelSize * nPixel / 2.0 + pixelSize / 2.0 + i * pixelSize);
//                Tbump->setY(-pixelSize * nPixel / 2.0 + pixelSize / 2.0 + j * pixelSize);
//                assemblyBumps->AddPlacedVolume(bumpLogicalVolume, *Tbump, Rbump);
//            }
//        }

//        G4ThreeVector *Tassembly = new G4ThreeVector();
//        assemblyBumps->MakeImprint(worldLV, *Tassembly, &fRotation);
//    }// end if (bumps == true)


//    G4double chipHeight = 300 * um;

//    if (chip == true) {
//        //Place an electronics chip behind the sensor

//        G4Box *chipBox = new G4Box("chip_box", 1.0 * pixelSize * nPixel / 2.0,
//                                   1.0 * pixelSize * nPixel / 2.0,
//                                   chipHeight / 2.0);
//        G4LogicalVolume *chipLogicalVolume = new G4LogicalVolume(chipBox,
//                G4Material::GetMaterial("G4_Si"),
//                "electronics_chip");
//        G4AssemblyVolume *chipassembly = new G4AssemblyVolume();

//        G4ThreeVector Tchip;

//        // rotate first
//        chipassembly->AddPlacedVolume(chipLogicalVolume, Tchip, &fRotation);

//        G4double ccChipSensor = sensorThickness /  2.0 + bumpHeight + chipHeight / 2.;
//        // translate relative to mother volume
//        Tchip.setZ(cos(fRotation.getTheta())*ccChipSensor);
//        Tchip.setY(-sin(fRotation.getTheta())*ccChipSensor);
//        // place "assembly"
////         G4cout << G4endl << "DEBUG rotation" << fRotation << G4endl << G4endl;
////         G4cout << G4endl << "DEBUG x angle" <<  fRotation.getPhi()*180/pi << G4endl << G4endl;
////         G4cout << G4endl << "DEBUG y angle" <<  fRotation.getPsi()*180/pi << G4endl << G4endl;
////         G4cout << G4endl << "DEBUG z angle" <<  fRotation.getTheta()*180/pi << " " << cos(pi) << " " << sin(pi) << G4endl << G4endl;
//        chipassembly->MakeImprint(worldLV, Tchip, new G4RotationMatrix());

//    }//end if (chip == true)

//    if (pcb == true) {
//        G4double pcbThickness = 0.5 * mm;
//        G4Box *pcbBox = new G4Box("pcb_box", 2.0 * pixelSize * nPixel / 2.0,
//                                  2.0 * pixelSize * nPixel / 2.0,
//                                  pcbThickness / 2.0);
//        G4LogicalVolume *pcbLogicalVolume = new G4LogicalVolume(pcbBox,
//                G4Material::GetMaterial("G4_SILICON_DIOXIDE"),
//                "printed_circuit_board");

//        G4AssemblyVolume *pcbAssembly = new G4AssemblyVolume();

//        G4ThreeVector Tpcb;

//        // rotate first
//        pcbAssembly->AddPlacedVolume(pcbLogicalVolume, Tpcb, &fRotation);

//        G4double ccPcbSensor = sensorThickness /  2.0 + bumpHeight + chipHeight + pcbThickness / 2.0;
//        // translate relative to mother volume
//        Tpcb.setZ(cos(fRotation.getTheta())*ccPcbSensor);
//        Tpcb.setY(-sin(fRotation.getTheta())*ccPcbSensor);
//        // place "assembly"
//        pcbAssembly->MakeImprint(worldLV, Tpcb, new G4RotationMatrix());


//        G4VisAttributes *pcb_visAtt
//            = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
//        pcb_visAtt->SetForceWireframe(true);
//        pcbLogicalVolume->SetVisAttributes(pcb_visAtt);

//    }//end if (pcb == true)


////-----------------------------------------------------------------Block
////--------------------------------------------------------------------------

//    if (block == true) {
//        G4Box *actBox = new G4Box("act_box", actX / 2.0,
//                                  actY / 2.0,
//                                  actZ / 2.0);
//        G4LogicalVolume *actLogicalVolume = new G4LogicalVolume(actBox,
//                G4Material::GetMaterial("G4_Fe"),
//                "activated_metal");
//        new G4PVPlacement(0,
//                          G4ThreeVector(0, 0, actPosZ - sensorThickness /  2.0 - actZ / 2.0),
//                          actLogicalVolume,
//                          "act_",
//                          worldLV,
//                          false,
//                          0);

//        G4VisAttributes *act_visAtt
//            = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
//        act_visAtt->SetForceWireframe(true);
//        actLogicalVolume->SetVisAttributes(act_visAtt);


//    }

////-------------------------------------------------------------------------------------------------------

//    //
//    // print parameters
//    //
//    G4cout << "\n------------------------------------------------------------\n" << G4endl;

//    G4cout << "Sensor: " << sensorThickness / um << "um of "
//           << sensorMaterial->GetName()
//           << " with " << nPixel << "x" << nPixel << " pixels at " << pixelSize / um << "um pitch" << G4endl;

//    if (filter == true)
//        G4cout << "Filter: " << filterThickness / um << "um of " << filterMaterial->GetName() << " at " << filterZ / um << "um " << G4endl;
//    else
//        G4cout << "Filter: disabled" << G4endl;

//    if (bumps == true)
//        G4cout << "Bumps: " << bumpMaterial->GetName() << " with r: " << bumpRadii / um << "um and h: " << bumpHeight / um << "um " << G4endl;
//    else
//        G4cout << "Bumps: disabled" << G4endl;




//    // create geometry here:
////    G4Material* fluorescenceMaterial1 = G4Material::GetMaterial("G4_Galactic");
////    G4Material* fluorescenceMaterial2 = G4Material::GetMaterial("G4_Galactic");

////    G4Box* fluorescenceS1           = new G4Box("flBox", 5*mm/2., 10*mm/2., 50*um/2);
////    G4LogicalVolume* fluorescenceLV1 = new G4LogicalVolume(fluorescenceS1,
////                                                         fluorescenceMaterial1,
////                                                         "flLV");


////    new G4PVPlacement(0,
////                  G4ThreeVector(5*mm/2,0,-0.53*cm),
////                  fluorescenceLV1,
////                  "fluorecence01_Cu",
////                  worldLV,
////                  false,
////                  0);


////     new G4PVPlacement(0,
////                   G4ThreeVector(5*mm/2,0,-0.53*cm),
////                   fluorescenceLV1,
////                   "fluorecence01_Cu",
////                   worldLV,
////                   false,
////                   0);



////    G4Box* fluorescenceS2           = new G4Box("flBox", 10*mm/2., 5*mm/2., 50*um/2);
////    G4LogicalVolume* fluorescenceLV2 = new G4LogicalVolume(fluorescenceS2,
////                                                         fluorescenceMaterial2,
////                                                         "flLV");

////    new G4PVPlacement(0,
////                  G4ThreeVector(0,5*mm/2,-0.5*cm),
////                  fluorescenceLV2,
////                  "fluorecence02_Ag",
////                  worldLV,
////                  false,
////                  0);



//    //
//    // Visualization attributes
//    //
//    G4VisAttributes *worldVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
//    worldVisAtt->SetVisibility(true);
//    //worldLV->SetVisAttributes(worldVisAtt);

//    worldLV->SetVisAttributes(G4VisAttributes::Invisible);

////     G4VisAttributes * flVisAtt1 = new G4VisAttributes(G4Colour(1.,0.,0.));
////     fluorescenceLV1->SetVisAttributes(flVisAtt1);

////     G4VisAttributes * flVisAtt2 = new G4VisAttributes(G4Colour(0.,1.,0.));
////     fluorescenceLV2->SetVisAttributes(flVisAtt2);

//    G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
//    simpleBoxVisAtt->SetVisibility(true);





//    //
//    // Always return the physical World
//    //
//    return worldPV;
//}
