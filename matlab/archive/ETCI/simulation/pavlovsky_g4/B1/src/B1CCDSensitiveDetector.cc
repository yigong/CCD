#include "B1CCDSensitiveDetector.hh"
#include "G4UnitsTable.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

B1CCDSensitiveDetector::B1CCDSensitiveDetector(G4String name)
:G4VSensitiveDetector(name)
{
}

B1CCDSensitiveDetector::~B1CCDSensitiveDetector()
{
}

void B1CCDSensitiveDetector::Initialize(G4HCofThisEvent *HCE)
{
}

G4bool B1CCDSensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
    G4double energy_dep = aStep->GetTotalEnergyDeposit();
    if(energy_dep==0.) return false;
    G4String name = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
    G4cout << "particle interaction: " << name << " energy of int: " << G4BestUnit(energy_dep, "energy");

}

void B1CCDSensitiveDetector::EndOfEvent(G4HCofThisEvent *)
{
}
