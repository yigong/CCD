#ifndef B1CCDHIT_HH
#define B1CCDHIT_HH

#include "G4VHit.hh"
#include "G4ThreeVector.hh"
#include "QFile"
#include "QTextStream"

class B1CCDHit : G4VHit
{
public:

    // post_step_volume_exit is 0 when the particle is in the ccd, 1 when it has exitted
    B1CCDHit(G4String particle_name, G4ThreeVector position, G4double kin_e,
                      G4double delta_e, G4double e_deposited, G4String post_step_process_name,
                      G4int generation_index, G4int parent_index, G4bool post_step_volume_exit, G4bool pre_step_volume_entry);
    ~B1CCDHit();
//    B1CCDHit(const B1CCDHit &right);
//    const B1CCDHit& operator=(const B1CCDHit &right);
//    G4int operator==(const B1CCDHit &right) const;

//    inline void *operator new(size_t);
//    inline void operator delete(void *aHit);

//    void Draw();
//    const std::map<G4String,G4AttDef>* GetAttDefs() const;
//    std::vector<G4AttValue>* CreateAttValues() const;
    void Print();
    void Print(QTextStream &out_stream);

private:
    // how many generations deep is the current track
    short int ge, pa; // to conserve memory space :(

    float edep_kev, ke_kev, de_kev; // in keV to conserve space in memory :(
    G4ThreeVector pos;
//    static std::map<G4String,G4AttDef> fAttDefs;
    G4String name, process;
    G4bool volexitcheck;
    G4bool volentcheck;

public:
    G4String GetParticleName() {return name;}
    float GetDE() { return de_kev;}
    float GetEdep() { return edep_kev;}
    G4ThreeVector GetPos() { return pos;}
    G4String GetParticleProcess() {return process;}
    float GetKE() { return ke_kev;}
    short int GetParentID() {return pa;}
    short int GetGenerationNumber() {return ge;}
    G4bool GetVolExitCheck() {return volexitcheck;}
    G4bool GetVolEntryCheck() { return volentcheck;}

};

typedef std::vector<B1CCDHit*> B1CCDHitCollection;



#endif // B1CCDHIT_HH
