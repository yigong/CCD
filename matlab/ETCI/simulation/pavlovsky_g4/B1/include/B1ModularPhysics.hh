#ifndef B1PHYSICSLIST_HH
#define B1PHYSICSLIST_HH

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class StepMax;

class B1ModularPhysics : public G4VModularPhysicsList
{
public:
    B1ModularPhysics();
    virtual ~B1ModularPhysics();

    void ConstructParticle();

    void SetCuts();
    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);

    void AddPhysicsList(const G4String& name);
    void ConstructProcess();

    void AddStepMax();
    StepMax* GetStepMaxProcess() {return stepMaxProcess;}

private:

  // these methods Construct physics processes and register them
//  void ConstructDecay();
//  void ConstructEM();

  G4VPhysicsConstructor* emPhysicsList;
  G4VPhysicsConstructor* decPhysicsList;
  G4VPhysicsConstructor*  particleList;
  StepMax* stepMaxProcess;

};

#endif // B1PHYSICSLIST_HH
