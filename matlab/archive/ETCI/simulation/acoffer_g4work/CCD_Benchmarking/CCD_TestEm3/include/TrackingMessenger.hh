/// \brief Definition of the TrackingMessenger class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TrackingMessenger_h
#define TrackingMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class TrackingAction;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackingMessenger: public G4UImessenger
{
  public:
    TrackingMessenger(TrackingAction*);
   ~TrackingMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    TrackingAction*   fTrackingAction;    
    G4UIcmdWithABool* fTrackingCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
