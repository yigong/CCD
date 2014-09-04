#include "B1CCDHit.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4AttValue.hh"
#include "G4AttDef.hh"
#include "G4AttCheck.hh"


B1CCDHit::B1CCDHit(G4String particle_name, G4ThreeVector position, G4double kin_e, G4double delta_e,
                   G4double e_deposited, G4String post_step_process_name, G4int generation_index, G4int parent_index,
                   G4bool post_step_volume_exit, G4bool pre_step_volume_entry)
{
    name = particle_name;
    de_kev = delta_e / keV;
    edep_kev = e_deposited / keV;
    pos = position;
    process = post_step_process_name;
    ke_kev = kin_e / keV;
    pa = parent_index;
    ge = (short int) generation_index;
    volexitcheck = (short int) post_step_volume_exit;
    volentcheck = (short int) pre_step_volume_entry;
}

B1CCDHit::~B1CCDHit()
{;}

//B1CCDHit::B1CCDHit(const B1CCDHit &right)
//  : G4VHit()
//{
//  edep = right.edep;
//  pos = right.pos;
//}

//const B1CCDHit& B1CCDHit::operator=(const B1CCDHit &right)
//{
//  edep = right.edep;
//  pos = right.pos;
//  return *this;
//}

//G4int B1CCDHit::operator==(const B1CCDHit &right) const
//{
//  return (this==&right) ? 1 : 0;
//}

//std::map<G4String,G4AttDef> B1CCDHit::fAttDefs;

//void B1CCDHit::Draw()
//{
//  // DO NOT IMPLEMENT - RP
//}

//const std::map<G4String,G4AttDef>* B1CCDHit::GetAttDefs() const
//{
//  // G4AttDefs have to have long life.  Use static member...
//  if (fAttDefs.empty()) {
//    fAttDefs["HitType"] =
//      G4AttDef("HitType","Type of hit","Physics","","G4String");
//  }
//  return &fAttDefs;
//}

//std::vector<G4AttValue>* B1CCDHit::CreateAttValues() const
//{
//  // Create expendable G4AttsValues for picking...
//  std::vector<G4AttValue>* attValues = new std::vector<G4AttValue>;
//  attValues->push_back
//    (G4AttValue("HitType","B1CCDHit",""));
//  return attValues;
//}

void B1CCDHit::Print()
{
    G4cout << name << "<" << pos.x() / um << ", " << pos.y() / um << ", " << pos.z() / um << "> "
           << " " << ke_kev << " " << de_kev << " " << edep_kev << " " << ge
           << " " << pa << " " << process << " " << " " << volentcheck <<  " "
           << volexitcheck << "\n";
}

void B1CCDHit::Print(QTextStream &out_stream)
{
    out_stream << name << "<" << pos.x() / um << ", " << pos.y() / um << ", " << pos.z() / um << "> "
           << " " << ke_kev << " " << de_kev << " " << edep_kev << " " << ge
           << " " << pa << " " << process << " " << " " << volentcheck <<  " "
           << volexitcheck << "\n";
}
