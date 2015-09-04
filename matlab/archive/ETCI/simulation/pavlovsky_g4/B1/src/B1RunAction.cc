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
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1EventAction.hh"
#include "B1SteppingAction.hh"
  // use of other actions 
  // - primary generator: to get info for printing about the primary
  // - event action: to get and reset accumulated energy sums
  // - stepping action: to get info about accounting volume 

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "QVector"
#include "QString"
#include "QTime"
#include "TH1D.h"
#include "TH2D.h"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "diffusion.hh"
#include "QTextStream"
#include "B1DetectorConstruction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    
  //initialize event cumulative quantities
  B1EventAction::Instance()->Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nofEvents = aRun->GetNumberOfEvent();
  if (nofEvents == 0) return;
  G4cout << "Number of events: " << nofEvents << "\n";
  // Compute dose
  //
  G4double energySum  = B1EventAction::Instance()->GetEnergySum();
  G4double energy2Sum = B1EventAction::Instance()->GetEnergy2Sum();
  G4double rms = energy2Sum - energySum*energySum/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

  // Run conditions
  //
  G4GeneralParticleSource* particleGun
    = B1PrimaryGeneratorAction::Instance()->GetParticleGun();
  G4String particleName
    = particleGun->GetParticleDefinition()->GetParticleName();
  G4double particleEnergy = particleGun->GetParticleEnergy();
  G4double mass = B1SteppingAction::Instance()->GetVolume()->GetMass();


  std::vector<B1CCDHitCollection*> *h_event_collection = B1EventAction::Instance()->GetEC();

  QVector<G4double>* energyArray = B1EventAction::Instance()->GetInteractionVals();
  G4cout << "Now";

  QString title = "Source Photons: "+QString::number(nofEvents,10);
  TH1D *hist = new TH1D("histogram",title.toAscii(),3000, 0., 1500);
  for(int i = 0; i < energyArray->size()-1; i++){
      G4double tmp = energyArray->at(i)/keV;
//      G4cout << tmp << "\n";
      hist->Fill( tmp );
  }
  QString out = QTime::currentTime().toString()+".root";
  hist->SaveAs(out.toAscii());
  delete(hist);

  if(h_event_collection->empty()) return;

  QFile *file = new QFile("tracks.out");
  file->open(QIODevice::WriteOnly);
  QTextStream out_stream(file);

  TH2D* loss_vs_dep = new TH2D("2dhistogram","Particle Energy Loss vs Energy Deposition", 1300, 0., 1300., 1300, 0., 1300.);
  TH2D* ccd_hist = new TH2D("histogram","CCD Diffused Tracks", 726*2, 0, 726-1, 1452*2, 0, 1454-1);
  diffusion* diff_track = new diffusion(ccd_hist,10.);

  QVector<G4String> process_list;
  QVector<int> process_counts;
  G4VPhysicalVolume* envLV = G4PhysicalVolumeStore::GetInstance()->GetVolume("Ccd");
  G4ThreeVector pos_0 = envLV->GetObjectTranslation();
  //ccd is offset from origin
  pos_0.setX( pos_0.x()-650./2. * um );

  for(int k=0; k < h_event_collection->size()-1; k++){
      G4cout << "Event #" << k;
      B1CCDHitCollection* h_hc = h_event_collection->at(k);

      if( h_hc->empty() ) continue;
      G4double init_e_loss = 0.;
      G4double ccd_e_dep = 0.;
      int j = 0;
      double energy=0;

      if( 0. == -h_hc->at(j)->GetDE() ) continue;
      short int parent_id = h_hc->at(0)->GetParentID();
      for(j=0; j < h_hc->size()-1; j++){
          G4ThreeVector pos = h_hc->at(j)->GetPos();
          energy=h_hc->at(j)->GetEdep();
          diff_track->diffuseInteraction((pos.x()-pos_0.x())/um,
                                         (pos.y()-pos_0.y())/um,
                                         (pos.z()-pos_0.z())/um,
                                         energy*1000. );
          h_hc->at(j)->Print(out_stream);
          bool is_direct_descendent = (h_hc->at(j)->GetParentID() == parent_id);
          if( is_direct_descendent ){
              init_e_loss +=-h_hc->at(j)->GetDE();
          }
          ccd_e_dep += h_hc->at(j)->GetEdep();

          int process_found = 0; int loc=0;
          G4String tmp_string = h_hc->at(j)->GetParticleProcess();
          for(int p=0; p<process_list.size()-1; p++){
              if( process_list.at(p) ==  tmp_string){ loc=p; process_found=1;}
          }
          if( !process_found ){
              process_list.push_back(tmp_string);
              process_counts.push_back(1);
          }else{
              process_counts.replace(loc, process_counts.at(loc)+1);
          }

      }
      loss_vs_dep->Fill(init_e_loss, ccd_e_dep );
          for(j=0; j < h_hc->size()-1; j++){
              if( init_e_loss < ccd_e_dep | init_e_loss > particleEnergy / keV ){
              h_hc->at(j)->Print(out_stream);
          }
      }
      out_stream << "\n";
  }
  for(int p=0; p<process_list.size()-1;p++){
      G4cout << process_list.at(p) << ": " << process_counts.at(p) << "\n";
  }
  file->close();
  ccd_hist->SaveAs("tracked.root");
  QString out_dep = "loss_vs_dep"+QTime::currentTime().toString()+".root";
  loss_vs_dep->SaveAs(out_dep.toAscii());
  delete loss_vs_dep, ccd_hist, diff_track, file;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
