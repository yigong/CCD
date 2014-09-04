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
/// \file electromagnetic/TestEm3/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "TrackingMessenger.hh"

#include "G4Track.hh"
#include "G4Positron.hh"
#include "G4PhysicalConstants.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det,RunAction* run,
                               EventAction* evt)
:fDetector(det), fRunAct(run), fEventAct(evt)
{
    fullChain = true;
    fTrackMessenger = new TrackingMessenger(this);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{
    delete fTrackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track )
{
  // Particle Tracking
    G4ParticleDefinition* particle = track->GetDefinition();
    G4String name   = particle->GetParticleName();
    fCharge = particle->GetPDGCharge();
    fMass   = particle->GetPDGMass();
    
    G4double Ekin = track->GetKineticEnergy();
    G4int ID      = track->GetTrackID();
    
    //count particles
    //
    fRunAct->ParticleCount(name, Ekin);
    
  // Energy flow initialisation for primary particle
  //
  if (track->GetTrackID() == 1) {
    G4int Idnow = 1;
    if (track->GetVolume() != fDetector->GetphysiWorld()) {
      // unique identificator of layer+absorber
      const G4VTouchable* touchable = track->GetTouchable();
      G4int absorNum = touchable->GetCopyNumber();
      G4int layerNum = touchable->GetReplicaNumber(1);
      Idnow = (fDetector->GetNbOfAbsor())*layerNum + absorNum;
    }
    
    G4double Eflow = track->GetKineticEnergy();
    if (track->GetDefinition() == G4Positron::Positron())
      Eflow += 2*electron_mass_c2; 
         
    //flux artefact, if primary vertex is inside the calorimeter   
    for (G4int pl=1; pl<=Idnow; pl++) {fRunAct->SumEnergyFlow(pl, Eflow);}
  } else {
    fRunAct->AddSecondaryTrack(track);
  }
  // Adding Radioactivity Ion Decay Tracking
    //fullChain: stop ion and print decay chain
    //
    if (fCharge > 2.) {
        G4Track* tr = (G4Track*) track;
        if (fullChain) tr->SetTrackStatus(fStopButAlive);
        if (ID == 1) fEventAct->AddDecayChain(name);
        else       fEventAct->AddDecayChain(" ---> " + name);
    }
  
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
    //keep only ions
    //
    if (fCharge < 3. ) return;
    
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();
    
    //get time
    //
    G4double time = track->GetGlobalTime();
    G4int ID = track->GetTrackID();
    if (ID == 1) fRunAct->PrimaryTiming(time);        //time of life of primary ion
    
    //energy and momentum balance (from secondaries)
    //
    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    size_t nbtrk = (*secondaries).size();
    if (nbtrk) {
        //there are secondaries --> it is a decay
        //
        //force 'single' decay
        if ((!fullChain)&&(ID > 1)) G4RunManager::GetRunManager()->AbortEvent();
        //
        //balance
        G4double EkinTot = 0.;
        G4ThreeVector Pbalance = - track->GetMomentum();
        for (size_t itr=0; itr<nbtrk; itr++) {
            G4Track* trk = (*secondaries)[itr];
            EkinTot += trk->GetKineticEnergy();
            //exclude gamma desexcitation from momentum balance
            if (trk->GetDefinition() != G4Gamma::Gamma())
                Pbalance += trk->GetMomentum();
        }
        G4double Pbal = Pbalance.mag();
        fRunAct->Balance(EkinTot,Pbal);
        analysis->FillH1(6,EkinTot);
        analysis->FillH1(7,Pbal);
    }
    
    //no secondaries --> end of chain
    //
    if (!nbtrk) {
        fRunAct->EventTiming(time);                     //total time of life
        analysis->FillH1(8,time);
    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

