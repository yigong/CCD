#include "B1PreDiffusedAnalysis.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "B1SteppingAction.hh"
#include "B1EventAction.hh"

#include "diffusion.hh"
#include "QTextStream"
#include "QVector"
#include "QString"
#include "QTime"
#include "TH1D.h"
#include "TH2D.h"


B1PreDiffusedAnalysis::B1PreDiffusedAnalysis(int numberofevents)
{
    n_of_events = numberofevents;
    if (n_of_events == 0) return;
    G4cout << "Number of events: " << n_of_events << "\n";

}

void B1PreDiffusedAnalysis::run(){


    // Compute dose
    //
    G4double energySum  = B1EventAction::Instance()->GetEnergySum();
    G4double energy2Sum = B1EventAction::Instance()->GetEnergy2Sum();
    G4double rms = energy2Sum - energySum*energySum/n_of_events;
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

    QString title = "Source Photons: "+QString::number(n_of_events,10);
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
