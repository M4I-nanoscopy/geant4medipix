#include "SteppingAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "G4RunManager.hh"


#include "G4SteppingManager.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

SteppingAction::SteppingAction(RunAction *run, EventAction *event) :
    G4UserSteppingAction(),
    fRunAction(run),
    fEventAction(event)
{
    G4RunManager *fRM = G4RunManager::GetRunManager();
    myDet = (DetectorConstructionBase *)(fRM->GetUserDetectorConstruction());
    thickness = myDet->GetSensorThickness();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
void SteppingAction::UserSteppingAction(const G4Step *aStep)
{
    // energy deposition per step
    G4double edep = aStep->GetTotalEnergyDeposit();
    if (edep <= 0.) return;

    //longitudinal profile of deposited energy
    G4ThreeVector prePoint  = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector postPoint = aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector point = prePoint + G4UniformRand() * (postPoint - prePoint);
    G4double r = point.z() + 0.5 * thickness; //FIXME get detector thickness
    

    // fill histograms
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    if( aStep->GetTrack()->GetVolume()->GetName() == "pixel_cell"){
      analysisManager->FillH1(1, r, edep);
      // total edep per event
      fEventAction->AddEdep(edep);
    }
}
