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
/// \file electromagnetic/TestEm3/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include <sstream>

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{ 
  fTestemDir = new G4UIdirectory("/testem/");
  fTestemDir->SetGuidance("UI commands specific to this example");
  
  fDetDir = new G4UIdirectory("/testem/det/");
  fDetDir->SetGuidance("detector construction commands");
  
  fSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeYZ",this);
  fSizeYZCmd->SetGuidance("Set tranverse size of the calorimeter");
  fSizeYZCmd->SetParameterName("Size",false);
  fSizeYZCmd->SetRange("Size>0.");
  fSizeYZCmd->SetUnitCategory("Length");
  fSizeYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fNbLayersCmd = new G4UIcmdWithAnInteger("/testem/det/setNbOfLayers",this);
  fNbLayersCmd->SetGuidance("Set number of layers.");
  fNbLayersCmd->SetParameterName("NbLayers",false);
  fNbLayersCmd->SetRange("NbLayers>0");
  fNbLayersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fNbAbsorCmd = new G4UIcmdWithAnInteger("/testem/det/setNbOfAbsor",this);
  fNbAbsorCmd->SetGuidance("Set number of Absorbers.");
  fNbAbsorCmd->SetParameterName("NbAbsor",false);
  fNbAbsorCmd->SetRange("NbAbsor>0");
  fNbAbsorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
  fAbsorCmd = new G4UIcommand("/testem/det/setAbsor",this);
  fAbsorCmd->SetGuidance("Set the absor nb, the material, the thickness.");
  fAbsorCmd->SetGuidance("  absor number : from 1 to NbOfAbsor");
  fAbsorCmd->SetGuidance("  material name");
  fAbsorCmd->SetGuidance("  thickness (with unit) : t>0.");
  //
  G4UIparameter* AbsNbPrm = new G4UIparameter("AbsorNb",'i',false);
  AbsNbPrm->SetGuidance("absor number : from 1 to NbOfAbsor");
  AbsNbPrm->SetParameterRange("AbsorNb>0");
  fAbsorCmd->SetParameter(AbsNbPrm);
  //
  G4UIparameter* MatPrm = new G4UIparameter("material",'s',false);
  MatPrm->SetGuidance("material name");
  fAbsorCmd->SetParameter(MatPrm);
  //    
  G4UIparameter* ThickPrm = new G4UIparameter("thickness",'d',false);
  ThickPrm->SetGuidance("thickness of absorber");
  ThickPrm->SetParameterRange("thickness>0.");
  fAbsorCmd->SetParameter(ThickPrm);
  //
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',false);
  unitPrm->SetGuidance("unit of thickness");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm->SetParameterCandidates(unitList);
  fAbsorCmd->SetParameter(unitPrm);
  //
  fAbsorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //
  // Create Source Positioning UI Command
  fSourcePositionCmd = new G4UIcommand("/testem/det/setSourcePosition",this);
  fSourcePositionCmd->SetGuidance("Set the x, y, z, and unit of the Emitting source Location");
  // Setup Parameters for Setting Source Volume Position
  G4UIparameter* xPrm = new G4UIparameter("x_position",'d',false);
  xPrm->SetGuidance("X Position of Source Volume, stay outside Cryostat Specs");
  fSourcePositionCmd->SetParameter(xPrm);
  
  G4UIparameter* yPrm = new G4UIparameter("y_position",'d',false);
  yPrm->SetGuidance("y Position of Source Volume, stay outside Cryostat Specs");
  fSourcePositionCmd->SetParameter(yPrm);
  
  G4UIparameter* zPrm = new G4UIparameter("z_position",'d',false);
  yPrm->SetGuidance("z Position of Source Volume, stay outside Cryostat Specs");
  fSourcePositionCmd->SetParameter(zPrm);
  
  G4UIparameter* unitPrm2 = new G4UIparameter("unit2",'s',false);
  unitPrm2->SetGuidance("unit of Source Position");
  G4String unitList2 = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm2->SetParameterCandidates(unitList2);
  fSourcePositionCmd->SetParameter(unitPrm2);
  
  fSourcePositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //
  /*fMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setField",this);
  fMagFieldCmd->SetGuidance("Define magnetic field.");
  fMagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  fMagFieldCmd->SetParameterName("Bz",false);
  fMagFieldCmd->SetUnitCategory("Magnetic flux density");
  fMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);*/
  
     
  fUpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  fUpdateCmd->SetGuidance("Update calorimeter geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fSizeYZCmd;
  delete fNbLayersCmd;
  delete fNbAbsorCmd;
  delete fAbsorCmd;
  delete fSourcePositionCmd;
  //delete fMagFieldCmd;
  delete fUpdateCmd;
  delete fDetDir;  
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fSizeYZCmd )
   { Detector->SetCalorSizeYZ(fSizeYZCmd->GetNewDoubleValue(newValue));}

  if( command == fNbLayersCmd )
   { Detector->SetNbOfLayers(fNbLayersCmd->GetNewIntValue(newValue));}

  if( command == fNbAbsorCmd )
   { Detector->SetNbOfAbsor(fNbAbsorCmd->GetNewIntValue(newValue));}
   
  if (command == fAbsorCmd)
   {
     G4int num; G4double tick;
     G4String unt, mat;
     std::istringstream is(newValue);
     is >> num >> mat >> tick >> unt;
     G4String material=mat;
     tick *= G4UIcommand::ValueOf(unt);
     Detector->SetAbsorMaterial (num,material);
     Detector->SetAbsorThickness(num,tick);
   }
  
    if (command == fSourcePositionCmd)
    {
      G4double x, y, z;
      //G4double tick2;
      G4String unt2;
      std::istringstream is2(newValue);
      is2 >> x >> y >> z >> unt2;
      x *= G4UIcommand::ValueOf(unt2);
      y *= G4UIcommand::ValueOf(unt2);
      z *= G4UIcommand::ValueOf(unt2);
      Detector->SetSourceVolumePosition(x, y, z);
      
    }
/*
  if( command == fMagFieldCmd )
   { Detector->SetMagField(fMagFieldCmd->GetNewDoubleValue(newValue));}
  */         
  if( command == fUpdateCmd )
   { Detector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
