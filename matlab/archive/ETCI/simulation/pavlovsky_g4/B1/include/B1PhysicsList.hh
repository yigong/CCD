#ifndef B1PHYSICSLIST_HH
#define B1PHYSICSLIST_HH

#include "G4VUserPhysicsList.hh"

class B1PhysicsList : public G4VUserPhysicsList
{
public:
    B1PhysicsList();
    virtual ~B1PhysicsList();

  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();

  void SetCuts();

private:

  // these methods Construct physics processes and register them
  void ConstructDecay();
  void ConstructEM();

};

#endif // B1PHYSICSLIST_HH
