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
// $Id$
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"

#include "B1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction* B1SteppingAction::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction* B1SteppingAction::Instance()
{
// Static acces function via G4RunManager 

  return fgInstance;
}      

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction()
: G4UserSteppingAction(),
  fVolume(0),
  fEnergy(0.)
{ 
  fgInstance = this;
  //Reset called before event action
//  fHitCollection = new B1CCDHitCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{ 
  fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
  if (volume != fVolume ) return;


  // collect energy and track length step by step
  G4double edep = step->GetTotalEnergyDeposit();
  fEnergy += edep;

  B1CCDHit *hit = new B1CCDHit(step->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName(),
                               step->GetTrack()->GetPosition(),
                               step->GetTrack()->GetKineticEnergy(),
                               step->GetDeltaEnergy(),
                               step->GetTotalEnergyDeposit(),
                               step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() ,
                               step->GetTrack()->GetTrackID(),
                               step->GetTrack()->GetParentID(),
                               fVolume == step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume(),
                               fVolume == volume);

  fHitCollection->push_back(hit);
//  G4cout << "Total step dep " <<step->GetTotalEnergyDeposit() / keV << "\n";
//  G4cout << "Particle Name " << step->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() << "\n";
//  G4cout << "Track id " << step->GetTrack()->GetTrackID() << "\n";
//  G4cout << "Parent ID " << step->GetTrack()->GetParentID() << "\n";
//  G4cout << "Kin E " << step->GetTrack()->GetKineticEnergy() / keV << "\n";
//  G4cout << "Is alive " << (step->GetTrack()->GetTrackStatus() == fAlive) << "\n";
//  G4cout << "Gen Proc " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << "\n";
//  G4cout << "Log Vol " << step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() <<"\n\n";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::Reset()
{
  fEnergy = 0.;

  //After reset, EventAction will own HC
  // DO NOT DELETE HERE!
  fHitCollection = new B1CCDHitCollection;

}

