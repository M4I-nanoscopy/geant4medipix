
#include "TrackingAction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

TrackingAction::TrackingAction()
    : G4UserTrackingAction()
{}

void TrackingAction::PostUserTrackingAction(const G4Track *aTrack)
{
    // count charge

//     G4TrackStatus status = aTrack->GetTrackStatus();

    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    //track length of primary particle or charged secondaries
    //
    G4double tracklen = 0.0;

    if ( aTrack->GetTrackID() == 1) {
//         G4cout << "DEBUG track length   " << (tracklen - 5000*um + 150*um)/um << " track id " << aTrack->GetTrackID() << G4endl;      
        tracklen = aTrack->GetTrackLength();
        tracklen -= 4000*um;
        analysisManager->FillH1(3,  tracklen);
    } else if (aTrack->GetParticleDefinition()->GetPDGCharge() != .0) {
    //else {
        tracklen = aTrack->GetTrackLength();
      //G4cout << "DEBUG track length   " << tracklen/um << " track id " << aTrack->GetTrackID() << G4endl;      
      analysisManager->FillH1(4, tracklen*um);
      }
    
}
