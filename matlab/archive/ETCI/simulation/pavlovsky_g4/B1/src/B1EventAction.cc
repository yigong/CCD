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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"

#include "B1RunAction.hh"
#include "B1SteppingAction.hh"
  // use of stepping action to get and reset accumulated energy  

#include "G4RunManager.hh"
#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction* B1EventAction::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction* B1EventAction::Instance()
{
// Static acces function via G4RunManager 

  return fgInstance;
}      

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction()
: G4UserEventAction(),
  fPrintModulo(1000),
  fEnergySum(0.),
  fEnergy2Sum(0.)
{ 
  fgInstance = this;
  fEnergyArray = new QVector<G4double>();
  fEventCollection = new std::vector<B1CCDHitCollection*>;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{ 
  fgInstance = 0;
  delete fEnergyArray;
  while( !fEventCollection->empty() ){
      B1CCDHitCollection* tmp = fEventCollection->back();
      fEventCollection->pop_back();
      while( !tmp->empty() ){
          B1CCDHit *tmp_hit = tmp->back();
          tmp->pop_back();
          delete tmp_hit;
      }
      delete tmp;
  }
  delete fEventCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event* event)
{  
  G4int eventNb = event->GetEventID();
  if (eventNb%fPrintModulo == 0) { 
    G4cout << "\n---> Begin of event: " << eventNb << G4endl;

    if( fEventCollection->size() ){
        //write to disk to avoid running out of ram
    }

  }
  // Reset accounted energy in stepping action
  B1SteppingAction::Instance()->Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event* event /*event*/)
{
  // accumulate statistics
  if((B1SteppingAction::Instance()->GetVolume()->GetName() == "Ccd") ){

      G4double energy = B1SteppingAction::Instance()->GetEnergy();
      if(energy!=0.0){
          B1CCDHitCollection *hc = B1SteppingAction::Instance()->GetHC();
          fEventCollection->push_back( hc );

          fEnergySum  += energy;
          fEnergy2Sum += energy*energy;
          fEnergyArray->push_back(energy);
      }else{
          B1CCDHitCollection* tmp = B1SteppingAction::Instance()->GetHC();
          while( !tmp->empty() ){
              B1CCDHit *tmp_hit = tmp->back();
              tmp->pop_back();
              delete tmp_hit;
          }
          delete tmp;
      }
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::Reset()
{
  //reset cumulative quantities
  //
  fEnergySum = 0.;
  fEnergy2Sum = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
