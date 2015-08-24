#ifndef B1CCDSENSITIVEDETECTOR_HH
#define B1CCDSENSITIVEDETECTOR_HH

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

class B1CCDSensitiveDetector : public G4VSensitiveDetector
{

public:
    B1CCDSensitiveDetector(G4String);
    ~B1CCDSensitiveDetector();

    void Initialize(G4HCofThisEvent*);
    G4bool  ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);
    void EndOfEvent(G4HCofThisEvent *);

private:
    // private vars
};

#endif // B1CCDSENSITIVEDETECTOR_HH
