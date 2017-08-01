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
/// \file src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
//
// $Id: PhysicsList.cc,v 1.37 2010-11-19 20:12:32 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////
//
// PhysicsList
//
// Created: 31.04.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
// 26.04.2007 Physics according to 8.3 Physics List (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
//

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

//radioactive decay, decay chain, de-exitation
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecay.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronElasticPhysicsHP.hh"
//#include "G4HadronElasticPhysicsLHEP.hh"

//#include "G4HadronQElasticPhysics.hh"

#include "G4ChargeExchangePhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4NeutronCrossSectionXS.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"


#include "G4IonPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmProcessOptions.hh"


#include "G4HadronInelasticQBBC.hh"

#include "FTFP_BERT_HP.hh"
#include "FTFP_BERT.hh"
#include "QGSP_FTFP_BERT.hh"

#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsFTF_BIC.hh"
#include "G4HadronPhysicsQGS_BIC.hh"
#include "G4HadronPhysicsQGSP_FTFP_BERT.hh"
#include "G4EmDNAPhysics.hh"

#include "G4IonPhysics.hh"

#include "G4LossTableManager.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4DNAGenericIonsManager.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
    G4LossTableManager::Instance();

    verboseLevel    = 0;

    fMessenger = new PhysicsListMessenger(this);

    G4LossTableManager::Instance();

    // Particles
    fParticleList = new G4DecayPhysics("decays", verboseLevel);
    // EM physics
    fEmPhysicsList = new G4EmStandardPhysics(verboseLevel);
    // radioactive decays
    fRadioactiveDecayList = new G4RadioactiveDecayPhysics(verboseLevel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PhysicsList::~PhysicsList()
{
    delete fMessenger;
    delete fParticleList;
    delete fEmPhysicsList;
    delete fRadioactiveDecayList;
    for (size_t i = 0; i < fHadronPhys.size(); i++) {
        delete fHadronPhys[i];
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::ConstructParticle()
{
    fParticleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::ConstructProcess()
{
    AddTransportation();
    fEmPhysicsList->ConstructProcess();
    fParticleList->ConstructProcess();

    // build radioactive decay
    fRadioactiveDecayList->ConstructProcess();
    G4RadioactiveDecay *radioactiveDecay = new G4RadioactiveDecay();
    radioactiveDecay->SetVerboseLevel(verboseLevel);
    radioactiveDecay->SetHLThreshold(nanosecond);
    radioactiveDecay->SetICM(true);
    radioactiveDecay->SetARM(false);

    for (size_t i = 0; i < fHadronPhys.size(); i++) {
        fHadronPhys[i]->ConstructProcess();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::AddPhysicsList(const G4String &name)
{
    if (verboseLevel > 0) {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
    }
    if (name == "emstandard_opt0") {

        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics();
    } else if (name == "emstandard_opt1") {

        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics_option1(verboseLevel);
    } else if (name == "emstandard_opt2") {

        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics_option2(verboseLevel);

    } else if (name == "emstandard_opt3") {

        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics_option3(verboseLevel);

    } else if (name == "livermore") {

        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmLivermorePhysics(verboseLevel);

    } else if (name == "penelope") {

//         delete fEmPhysicsList;
        fEmPhysicsList = new G4EmPenelopePhysics(verboseLevel);

    } else if (name == "FTFP_BERT_EMV") {

        AddPhysicsList("emstandard_opt1");
        AddPhysicsList("FTFP_BERT");

    } else if (name == "FTFP_BERT_EMX") {

        AddPhysicsList("emstandard_opt2");
        AddPhysicsList("FTFP_BERT");

    } else if (name == "FTFP_BERT") {

        SetBuilderList1();
        fHadronPhys.push_back(new G4HadronPhysicsFTFP_BERT());

    } else if (name == "FTF_BIC") {

        SetBuilderList0();
        fHadronPhys.push_back(new G4HadronPhysicsFTF_BIC());
        fHadronPhys.push_back(new G4NeutronCrossSectionXS(verboseLevel));

    } else if (name == "QBBC") {

        AddPhysicsList("emstandard_opt2");
        SetBuilderList3();
        fHadronPhys.push_back(new G4HadronInelasticQBBC());

    } else if (name == "QGSP_BERT") {

        SetBuilderList1();
        fHadronPhys.push_back(new G4HadronPhysicsQGSP_BERT());

    } else if (name == "QGSP_FTFP_BERT") {

        SetBuilderList1();
        fHadronPhys.push_back(new G4HadronPhysicsQGSP_FTFP_BERT());

    } else if (name == "QGSP_BERT_EMV") {

        AddPhysicsList("emstandard_opt1");
        AddPhysicsList("QGSP_BERT");

    } else if (name == "QGSP_BERT_EMX") {

        AddPhysicsList("emstandard_opt2");
        AddPhysicsList("QGSP_BERT");

    } else if (name == "QGSP_BERT_HP") {

        SetBuilderList1(true);
        fHadronPhys.push_back(new G4HadronPhysicsQGSP_BERT_HP());

    } else if (name == "QGSP_BIC") {

        SetBuilderList0();
        fHadronPhys.push_back(new G4HadronPhysicsQGSP_BIC());

    } else if (name == "QGSP_BIC_EMY") {

        AddPhysicsList("emstandard_opt3");
        SetBuilderList0();
        fHadronPhys.push_back(new G4HadronPhysicsQGSP_BIC());

    } else if (name == "QGS_BIC") {

        SetBuilderList0();
        fHadronPhys.push_back(new G4HadronPhysicsQGS_BIC());
        fHadronPhys.push_back(new G4NeutronCrossSectionXS(verboseLevel));

    } else if (name == "QGSP_BIC_HP") {

        SetBuilderList0(true);
        fHadronPhys.push_back(new G4HadronPhysicsQGSP_BIC_HP());

    } else {

        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
               << " is not defined"
               << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList0(G4bool flagHP)
{
    fHadronPhys.push_back(new G4EmExtraPhysics(verboseLevel));
    if (flagHP) {
        fHadronPhys.push_back(new G4HadronElasticPhysicsHP(verboseLevel));
    } else {
        fHadronPhys.push_back(new G4HadronElasticPhysics(verboseLevel));
    }
    fHadronPhys.push_back(new G4StoppingPhysics(verboseLevel));
    fHadronPhys.push_back(new G4IonBinaryCascadePhysics(verboseLevel));
    fHadronPhys.push_back(new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList1(G4bool flagHP)
{
    fHadronPhys.push_back(new G4EmExtraPhysics(verboseLevel));
    if (flagHP) {
        fHadronPhys.push_back(new G4HadronElasticPhysicsHP(verboseLevel));
    } else {
        fHadronPhys.push_back(new G4HadronElasticPhysics(verboseLevel));
    }
    fHadronPhys.push_back(new G4StoppingPhysics(verboseLevel));
    fHadronPhys.push_back(new G4IonPhysics(verboseLevel));
    fHadronPhys.push_back(new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList2(G4bool addStopping)
{

    fHadronPhys.push_back(new G4EmExtraPhysics(verboseLevel));
    fHadronPhys.push_back(new G4HadronElasticPhysics(verboseLevel));
    if (addStopping) {
        fHadronPhys.push_back(new G4StoppingPhysics(verboseLevel));
    }
    fHadronPhys.push_back(new G4IonPhysics(verboseLevel));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList3()
{
    fHadronPhys.push_back(new G4EmExtraPhysics(verboseLevel));
    RegisterPhysics(new G4HadronElasticPhysicsXS(verboseLevel));
    fHadronPhys.push_back(new G4StoppingPhysics(verboseLevel));
    fHadronPhys.push_back(new G4IonBinaryCascadePhysics(verboseLevel));
    fHadronPhys.push_back(new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysicsList::SetBuilderList4()
{
    fHadronPhys.push_back(new G4EmExtraPhysics(verboseLevel));
    fHadronPhys.push_back(new G4HadronElasticPhysics(verboseLevel));
    fHadronPhys.push_back(new G4StoppingPhysics(verboseLevel));
    fHadronPhys.push_back(new G4IonPhysics(verboseLevel));
    fHadronPhys.push_back(new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
void PhysicsList::List()
{
    G4cout << "### PhysicsLists available: FTFP_BERT FTFP_BERT_EMV FTFP_BERT_EMX FTF_BIC"
           << G4endl;
    G4cout << "                            LHEP LHEP_EMV QBBC QGS_BIC QGSP"
           << G4endl;
    G4cout << "                            QGSC_BERT QGSP_BERT QGSP_BERT_EMV QGSP_BIC_EMY"
           << G4endl;
    G4cout << "                            QGSP_BERT_EMX QGSP_BERT_HP QGSP_BIC QGSP_BIC_HP"
           << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

