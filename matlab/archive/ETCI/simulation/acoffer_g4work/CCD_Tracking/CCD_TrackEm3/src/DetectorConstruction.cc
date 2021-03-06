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
/// \file electromagnetic/TestEm3/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
//#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RegionStore.hh"

#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:fDefaultMaterial(0),  fSolidWorld(0),    fLogicWorld(0), fPhysiWorld(0),
 sourceContainerMaterial(0),fSourceVolume_x(-10.0*cm),fSourceVolume_y(0),fSourceVolume_z(0),
 fSourceVolumeTubs(0), fLogicSource(0), fPhysiSource(0),
 cryoMaterial(0),     cryoBox(0),        cryoLog(0),     cryoPhys(0),
 cryoCutoutMaterial(0), cryoCutoutBox(0),	cryoCutoutLog(0),  	cryoCutoutPhys(0),
 frontCutoutTubs(0), cryoBoxSub_1(0), cryoBoxSub(0), 
 connectorBox(0),	connectorLog(0),  	connectorBoxPhys(0),
 connectorCutoutBox(0),	connectorCutoutLog(0),  connectorCutoutPhys(0),
 ifBox(0), 		ifLog(0),  		ifBoxPhys(0),
 ifCutoutBox(0),	ifCutoutLog(0),  	ifCutoutPhys(0),
 fSolidCalor(0),fLogicCalor(0),fPhysiCalor(0),
 fSolidLayer(0),fLogicLayer(0),fPhysiLayer(0),
 CryoWallAbsorberStandOff(0), CalorOffset()
 
{
  // default parameter values of the calorimeter
  fNbOfAbsor = 3;
  fAbsorThickness[1] = 2.0*um;    // Pixel Face
  fAbsorThickness[2] = 641.0*um;  // CCD Depeated Region
  fAbsorThickness[3] = 7.0*um;    //
  
  // Region Names
  //constructed = false;
  absorName[0] = "Layer-A";
  absorName[1] = "Layer-B";
  absorName[2] = "Layer-C";
  
  //calorRegionName = "Layer-A";
  
  // CryoStat Dimentions of Use:
  
    /*
  fAbsorThickness[1] = 1.3*mm;
  fAbsorThickness[2] = 2.0*cm;
  fAbsorThickness[3] = 650*um;
  fAbsorThickness[4] = 2.0*cm;
  fAbsorThickness[5] = 5.0*mm;
   */
  fNbOfLayers        = 1;
  fCalorSizeYZ       = 3.*cm;
  
  // Source Volume Info
  //fSourceVolume_x         = -15.0*cm;
  //fSourceVolume_y         = 0.0*cm;
  //fSourceVolume_z         = 0.0*cm;
  fSourceVolumeDiameter   = 2.54*cm;
  fSourceVolumeThickness  = 6.35*mm;
  fSourceWindowThickness  = 2.77*mm;
  
  // Definition of cryostat external dims
  cryo_hx           = 6.1*cm;
  cryo_hy           = 17.3*cm;
  cryo_hz           = 23.4*cm;
  
  cryo_thickness_x  = 0.5*cm;
  cryo_thickness_y  = 0.7*cm;
  cryo_thickness_z  = cryo_thickness_y;
  
  ComputeCalorParameters();
  fCalorOffset = -fCalorThickness/2;
  
    // materials
  DefineMaterials();
  SetWorldMaterial("Air20");  // defines fDefaultMaterial (World and Non-Vac Volumes)
  SetSourceMaterial("HDPE");
  
  SetAbsorMaterial(1,"Silicon");
  SetAbsorMaterial(2,"Silicon" );
  SetAbsorMaterial(3,"Silicon" );
  
  SetCryoMaterial("Aluminium");
  SetCryoCutoutMaterial("Vac");
  
  /*
  SetAbsorMaterial(2,"Beam");
  SetAbsorMaterial(3,"Silicon");
  SetAbsorMaterial(4,"Beam");
  SetAbsorMaterial(5,"Aluminium");
  */
  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // This function illustrates the possible ways to define materials using 
  // G4 database on G4Elements
  G4NistManager* manager = G4NistManager::Instance();
  manager->SetVerbose(0);
  //
  // define Elements
  //
  
  G4double z,a;
  G4Element* H  = manager->FindOrBuildElement(1);
  G4Element* C  = manager->FindOrBuildElement(6);
  G4Element* N  = manager->FindOrBuildElement(7);
  G4Element* O  = manager->FindOrBuildElement(8);
  G4Element* Si = manager->FindOrBuildElement(14);
  G4Element* Ge = manager->FindOrBuildElement(32);
  G4Element* Sb = manager->FindOrBuildElement(51);
  G4Element* I  = manager->FindOrBuildElement(53);
  G4Element* Cs = manager->FindOrBuildElement(55);
  G4Element* Pb = manager->FindOrBuildElement(82);
  G4Element* Bi = manager->FindOrBuildElement(83);

  //
  // define an Element from isotopes, by relative abundance
  //
  //G4int iz, n;                       //iz=number of protons  in an isotope;
                                     // n=number of nucleons in an isotope;
  G4int   ncomponents;                                     
  //G4double abundance;
  /*
  G4Isotope* U5 = new G4Isotope("U235", iz=92, n=235, a=235.01*g/mole);
  G4Isotope* U8 = new G4Isotope("U238", iz=92, n=238, a=238.03*g/mole);

  G4Element* U  = new G4Element("enriched Uranium", "U", ncomponents=2);
  U->AddIsotope(U5, abundance= 90.*perCent);
  U->AddIsotope(U8, abundance= 10.*perCent);
  */

  //
  // define simple materials
  //
  G4double density;

  new G4Material("liquidH2",    z=1.,  a= 1.008*g/mole,  density= 70.8*mg/cm3);
  new G4Material("Aluminium",   z=13., a= 26.98*g/mole,  density= 2.700*g/cm3);
  new G4Material("Silicon",     z=14., a=28.099*g/mole,  density=2.33*g/cm3);
  new G4Material("Titanium",    z=22., a= 47.867*g/mole, density= 4.54*g/cm3);
  new G4Material("Iron",        z=26., a= 55.85*g/mole,  density= 7.870*g/cm3);
  new G4Material("Copper",      z=29., a= 63.55*g/mole,  density= 8.960*g/cm3);
  new G4Material("Tungsten",    z=74., a= 183.85*g/mole, density= 19.30*g/cm3);
  new G4Material("Gold",        z=79., a= 196.97*g/mole, density= 19.32*g/cm3);
  new G4Material("Uranium",     z=92., a= 238.03*g/mole, density= 18.95*g/cm3);

  //
  // define a material from elements.   case 1: chemical molecule
  //
  G4int natoms;

  G4Material* H2O = 
  new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  H2O->SetChemicalFormula("H_2O");
  
  G4Material* CH = 
  new G4Material("Polystyrene", density= 1.032*g/cm3, ncomponents=2);
  CH->AddElement(C, natoms=1);
  CH->AddElement(H, natoms=1);
  
  G4Material* CH2 =
  new G4Material("HDPE", density= 0.952*g/cm3, ncomponents=2);
  CH2->AddElement(C, natoms=1);
  CH2->AddElement(H, natoms=2);

  G4Material* Sci = 
  new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Sci->AddElement(C, natoms=9);
  Sci->AddElement(H, natoms=10);
  
  Sci->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  G4Material* Lct =
  new G4Material("Lucite", density= 1.185*g/cm3, ncomponents=3);
  Lct->AddElement(C, 59.97*perCent);
  Lct->AddElement(H, 8.07*perCent);
  Lct->AddElement(O, 31.96*perCent);

  G4Material* Sili = 
  new G4Material("Silicon", density= 2.330*g/cm3, ncomponents=1);
  Sili->AddElement(Si, natoms=1);

  G4Material* SiO2 = 
  new G4Material("quartz", density= 2.200*g/cm3, ncomponents=2);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O , natoms=2);

  G4Material* G10 = 
  new G4Material("NemaG10", density= 1.700*g/cm3, ncomponents=4);
  G10->AddElement(Si, natoms=1);
  G10->AddElement(O , natoms=2);
  G10->AddElement(C , natoms=3);
  G10->AddElement(H , natoms=3);

  G4Material* CsI = 
  new G4Material("CsI", density= 4.534*g/cm3, ncomponents=2);
  CsI->AddElement(Cs, natoms=1);
  CsI->AddElement(I , natoms=1);
  CsI->GetIonisation()->SetMeanExcitationEnergy(553.1*eV);

  G4Material* BGO = 
  new G4Material("BGO", density= 7.10*g/cm3, ncomponents=3);
  BGO->AddElement(O , natoms=12);
  BGO->AddElement(Ge, natoms= 3);
  BGO->AddElement(Bi, natoms= 4);

  //SiNx
  density= 3.1 *g/cm3;
  G4Material* SiNx= new G4Material("SiNx", density, ncomponents=3);
  SiNx-> AddElement(Si, 300);
  SiNx-> AddElement(N, 310);
  SiNx-> AddElement(H, 6);
  
  //
  // define gaseous materials using G4 NIST database 
  //
  G4double fractionmass;
  // Room Temp Air
  G4Material* Air = manager->FindOrBuildMaterial("G4_AIR");
  manager->ConstructNewGasMaterial("Air20","G4_AIR",293.*kelvin,1.*atmosphere);
  // Vacume Volume
  G4double detTemp = 130.0*kelvin; 
  G4Material* Vacuum = new G4Material("Vac", 14., 28.*g/mole, 1.e-7*g/cm3,
                                      kStateGas, detTemp, 1.e-6*pascal);
  
  //
  // define a material from elements and others materials (mixture of mixtures)
  //

  G4Material* Lead = new G4Material("Lead", density= 11.35*g/cm3, ncomponents=1);
  Lead->AddElement(Pb, fractionmass=1.0);

  G4Material* LeadSb = new G4Material("LeadSb", density= 11.35*g/cm3, ncomponents=2);
  LeadSb->AddElement(Sb, fractionmass=4.*perCent);
  LeadSb->AddElement(Pb, fractionmass=96.*perCent);

  G4Material* Aerog = new G4Material("Aerogel", density= 0.200*g/cm3, ncomponents=3);
  Aerog->AddMaterial(SiO2, fractionmass=62.5*perCent);
  Aerog->AddMaterial(H2O , fractionmass=37.4*perCent);
  Aerog->AddElement (C   , fractionmass= 0.1*perCent);

  //
  // examples of gas in non STP conditions
  //
  G4double temperature, pressure;
  
  G4Material* CO2 = 
  new G4Material("CarbonicGas", density= 27.*mg/cm3, ncomponents=2,
                 kStateGas, temperature= 325.*kelvin, pressure= 50.*atmosphere);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);

  G4Material* steam = 
  new G4Material("WaterSteam", density= 1.0*mg/cm3, ncomponents=1,
                  kStateGas, temperature= 273*kelvin, pressure= 1*atmosphere);
  steam->AddMaterial(H2O, fractionmass=1.);
  
  new G4Material("ArgonGas", z=18, a=39.948*g/mole, density= 1.782*mg/cm3,
                 kStateGas, 273.15*kelvin, 1*atmosphere);
  //
  // examples of vacuum
  //

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1., a=1.008*g/mole, density,
                             kStateGas,temperature,pressure);

  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  //temperature = STP_Temperature;         //from PhysicalConstants.h
  G4Material* beam = 
      new G4Material("Beam", density, ncomponents=1,
                         kStateGas,detTemp,pressure);
      beam->AddMaterial(Air, fractionmass=1.);

  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
  fLayerThickness = 0.;
  for (G4int iAbs=1; iAbs<=fNbOfAbsor; iAbs++) {
    fLayerThickness += fAbsorThickness[iAbs];
  }
  fCalorThickness = fNbOfLayers*fLayerThickness;     
    //fWorldSizeX = 2.2*fCalorThickness;
    //fWorldSizeYZ = 2.2*fCalorSizeYZ;
  fCalorOffset = fCalorThickness/2;
  CryoWallAbsorberStandOff = (fCalorThickness/2 +   25.675*mm);
  fWorldSizeX = 2.2*fCalorThickness + 2.2*cryo_hx + 2.2*std::abs(fSourceVolume_x);
  fWorldSizeYZ = 2.2*fCalorSizeYZ + 2.2*cryo_hz + 2.2*std::abs(fSourceVolume_y)+ 2.2*std::abs(fSourceVolume_z);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{
  // complete the Calor parameters definition
  ComputeCalorParameters();

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  //G4RegionStore::GetInstance()->clear();
  
  //
  // World
  //

  fSolidWorld = new G4Box("World",                                //its name
                   fWorldSizeX/2,fWorldSizeYZ/2,fWorldSizeYZ/2);        //its size

  fLogicWorld = new G4LogicalVolume(fSolidWorld,                //its solid
                                   fDefaultMaterial,        //its material
                                   "World");                //its name

  fPhysiWorld = new G4PVPlacement(0,                        //no rotation
                                   G4ThreeVector(),        //at (0,0,0)
                                 fLogicWorld,                //its fLogical volume
                                 "World",                //its name
                                 0,                        //its mother  volume
                                 false,                        //no boolean operation
                                 0);                        //copy number
//
  // Source Volume
//
  fSource_inner_rad = 0.*mm;
   fSourceVolumeTubs = new G4Tubs("sourceVolume",
                                 fSource_inner_rad,
                                 fSourceVolumeDiameter/4.0,
                                 fSourceVolumeThickness/2.0,
                                 0.*deg,
                                 360.*deg);
  fLogicSource = new G4LogicalVolume(fSourceVolumeTubs,               // its solid
                                     sourceContainerMaterial,    // its material
                                     "sourceVolume");            // its name
  G4RotationMatrix* yRotSourceVolume = new G4RotationMatrix();
  yRotSourceVolume->rotateY(M_PI/2.0*rad);
  fPhysiSource = new G4PVPlacement(yRotSourceVolume,                              // rotation,
                                  G4ThreeVector(fSourceVolume_x, //-fSourceVolumeThickness/2+fSourceWindowThickness),
                                                fSourceVolume_y,
                                                fSourceVolume_z), // its placement
                                  fLogicSource,                   // its fLogical volume
                                  "sourceVolume",                 // its name
                                  fLogicWorld,                        // its mother volume
                                  false,                          // no boolean operation
                                  0);                             // copy number 
                                  
  
                                     
  
//
// NEW Cryostat Box - 'Realistic'
//
    /// ..............................
    /// Start creating geometrical volumes
    /// ..............................
   // X is depth, Y is height and Z is width
   // Can be used to Shift the center of the Cryostat off of 0,0,0
    world_hx  = 0.50*m;
    world_hy  = 0.50*m;
    world_hz  = 0.50*m;
  
    /// Start Cryostat definitions
    
    cryoBox = new G4Box("Cryo", cryo_hx/2.0, cryo_hy/2.0, cryo_hz/2.0);
    
    // Definition of vacuum volume and Cryo-Cutout Volume 
    cryo_cutout_hx = cryo_hx-2*cryo_thickness_x;
    cryo_cutout_hy = cryo_hy-2*cryo_thickness_y;
    cryo_cutout_hz = cryo_hz-2*cryo_thickness_z;
    
    cryoCutoutBox = new G4Box("CryoCutout", cryo_cutout_hx/2.0, cryo_cutout_hy/2.0, cryo_cutout_hz/2.0);
    
    // Definition of front cutout
    front_cutout_inner_rad = 0.*cm;
    front_cutout_outer_rad = 4.99*cm;
    front_cutout_hx = 0.4*cm;
    
    frontCutoutTubs = new G4Tubs("frontCutout",
                                         front_cutout_inner_rad,
                                         front_cutout_outer_rad,
                                         front_cutout_hx/2.0,
                                         0.*deg,
                                         360.*deg);
    
    G4ThreeVector xTransFrontCutout_1(-0.5*(cryo_hx-front_cutout_hx)-0.01*cm, 0.*cm, 0.*cm);
    G4RotationMatrix* yRotFrontCutout = new G4RotationMatrix();
    yRotFrontCutout->rotateY(M_PI/2.0*rad);
    
    cryoBoxSub_1 = new G4SubtractionSolid("CyroCutoutComplete_1",
                                                              cryoBox,
                                                              frontCutoutTubs,
                                                              yRotFrontCutout,
                                                              xTransFrontCutout_1);
    
    G4ThreeVector xTransFrontCutout(0.5*(cryo_hx-front_cutout_hx)+0.01*cm, 0.*cm, 0.*cm);
    cryoBoxSub = new G4SubtractionSolid("CyroCutoutComplete",
                                                            cryoBoxSub_1,
                                                            frontCutoutTubs,
                                                            yRotFrontCutout,
                                                            xTransFrontCutout);
    /// Start electronics boxes defintion
    
    // Define connector box
    connector_hx  = 5.1*cm;
    connector_hy  = 16.3*cm;
    connector_hz  = 2.526*cm;
    
    connectorBox = new G4Box("ConnectorBox", connector_hx/2.0, connector_hy/2.0, connector_hz/2.0);
    
    // Define connector box cutout
    connector_thickness = 0.28*cm;
    connector_cutout_hx  = connector_hx-2.0*connector_thickness;
    connector_cutout_hy  = connector_hy-2.0*connector_thickness;
    connector_cutout_hz  = connector_hz-2.0*connector_thickness;
    
    connectorCutoutBox = new G4Box("ConnectorCutoutBox",
                                   connector_cutout_hx/2.0,
                                   connector_cutout_hy/2.0,
                                   connector_cutout_hz/2.0);
    
    // Define interface board box
    G4double if_hx  = 8.3*cm;
    G4double if_hy  = 19.5*cm;
    G4double if_hz  = 6.0*cm;
    
    ifBox = new G4Box("IFBox", if_hx/2.0, if_hy/2.0, if_hz/2.0);
    
    // Define interface board box
    G4double if_thickness = 0.28*cm;
    G4double if_cutout_hx  = if_hx - 2.0*if_thickness;
    G4double if_cutout_hy  = if_hy - 2.0*if_thickness;
    G4double if_cutout_hz  = if_hz - 2.0*if_thickness;
    
    ifCutoutBox = new G4Box("IFCutoutBox",
                            if_cutout_hx/2.0,
                            if_cutout_hy/2.0,
                            if_cutout_hz/2.0);
  
    /// ..............................
    /// Start Creating Logical Volumes
    /// ..............................
    //worldLog = new G4LogicalVolume(worldBox, Air, "World");
    cryoLog = new G4LogicalVolume(cryoBoxSub, cryoMaterial, "Cryo");
    cryoCutoutLog = new G4LogicalVolume(cryoCutoutBox, cryoCutoutMaterial, "CryoCutout");
    connectorLog = new G4LogicalVolume(connectorBox, cryoMaterial, "ConnectorBox");
    connectorCutoutLog = new G4LogicalVolume(connectorCutoutBox, fDefaultMaterial, "ConnectorBoxInternal");
    ifLog = new G4LogicalVolume(ifBox, cryoMaterial, "IFBox");
    ifCutoutLog = new G4LogicalVolume(ifCutoutBox, fDefaultMaterial, "IFBoxInternal");
    
    /// ..............................
    /// Start creating phsyical volumes
    /// ..............................
    G4double pos_x = 0.*world_hx;
    G4double pos_y = 0.*world_hy;
    G4double pos_z = 0.*world_hz;
    /*
    physWorld = new G4PVPlacement(0,
                                  G4ThreeVector(),       //at (0,0,0)
                                  fLogicWorld,            //its logical volume
                                  "World",               //its name
                                  0,                     //its mother  volume
                                  false,                 //no boolean operation
                                  0);                    //copy number
    */
    cryoPhys = new G4PVPlacement(0, //no rotation
                                 G4ThreeVector(pos_x,pos_y,pos_z),
                                 cryoLog,
                                 "Cryo",
                                 fLogicWorld,
                                 false,
                                 0);
    
    cryoCutoutPhys = new G4PVPlacement(0,
                                       G4ThreeVector(0,0,0),
                                       cryoCutoutLog,
                                       "CryoCutout",
                                       cryoLog,
                                       false,
                                       0);
    
    connectorBoxPhys = new G4PVPlacement(0, //no rotation
                                         G4ThreeVector(pos_x,pos_y,pos_z+0.5*(cryo_hz+connector_hz)),
                                         connectorLog,
                                         "ConnectorBox",
                                         fLogicWorld,
                                         false,
                                         0);
    
    
    connectorCutoutPhys = new G4PVPlacement(0, //no rotation
                                            G4ThreeVector(0.,0.,0.),
                                            connectorCutoutLog,
                                            "ConnectorCutout",
                                            connectorLog,
                                            false,
                                            0);
    
    ifBoxPhys = new G4PVPlacement(0, //no rotation
                                  G4ThreeVector(pos_x,
                                                pos_y-0.5*(if_hy-connector_hy)+0.28*cm,
                                                pos_z+0.5*(cryo_hz+if_hz)+connector_hz),
                                  ifLog,
                                  "ConnectorBox",
                                  fLogicWorld,
                                  false,
                                  0);
    
    
    ifCutoutPhys = new G4PVPlacement(0, //no rotation
                                     G4ThreeVector(0.,0.,0.),
                                     ifCutoutLog,
                                     "ConnectorCutout",
                                     ifLog,
                                     false,
                                     0);

    
    
    
  //
  // Calorimeter
  //

  fSolidCalor = new G4Box("Calorimeter",                                     //its name
                           fCalorThickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2);//size

  fLogicCalor = new G4LogicalVolume(fSolidCalor,                //its solid
                                         fDefaultMaterial,        //its material
                                         "Calorimeter");        //its name

  fPhysiCalor = new G4PVPlacement(0,                        //no rotation
                                 G4ThreeVector(fCalorOffset,0,0),        //at (0,0,0)
                                 fLogicCalor,                //its fLogical volume
                                 "Calorimeter",                //its name
                                 cryoCutoutLog,                //its mother  volume
                                 false,                        //no boolean operation
                                 0);                        //copy number

  //
  // Layers
  //

  fSolidLayer = new G4Box("Layer",                                      //its name
                       fLayerThickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); //size

  fLogicLayer = new G4LogicalVolume(fSolidLayer,                //its solid
                                   fDefaultMaterial,        //its material
                                   "Layer");                //its name
  if (fNbOfLayers > 1)
    fPhysiLayer = new G4PVReplica("Layer",                //its name
                                       fLogicLayer,                //its fLogical volume
                                       fLogicCalor,                //its mother
                                 kXAxis,                //axis of replication
                                 fNbOfLayers,                //number of replica
                                 fLayerThickness);        //witdth of replica
  else
    fPhysiLayer = new G4PVPlacement(0,                        //no rotation
                                   G4ThreeVector(),        //at (0,0,0)
                                   fLogicLayer,                //its fLogical volume
                                   "Layer",                //its name
                                   fLogicCalor,                //its mother  volume
                                   false,                //no boolean operation
                                   0);                        //copy number

  //
  // Absorbers
  //

  G4double xfront = -0.5*fLayerThickness;
  for (G4int k=1; k<=fNbOfAbsor; k++) {
    fSolidAbsor[k] = new G4Box("Absorber",                //its name
                              fAbsorThickness[k]/2,fCalorSizeYZ/2,fCalorSizeYZ/2);

    fLogicAbsor[k] = new G4LogicalVolume(fSolidAbsor[k],    //its solid
                                        fAbsorMaterial[k], //its material
                                        fAbsorMaterial[k]->GetName());

    G4double xcenter = xfront+0.5*fAbsorThickness[k];
    xfront += fAbsorThickness[k];
    fPhysiAbsor[k] = new G4PVPlacement(0,                   //no rotation
                         G4ThreeVector(xcenter,0.,0.),      //its position
                         fLogicAbsor[k],                    //its logical volume        
                         fAbsorMaterial[k]->GetName(),      //its name
                         fLogicLayer,                       //its mother
                         false,                             //no boulean operat
                         k);                                //copy number

  }
  
  //
  // Regions
  //
  /*
  // Absorber as Regions (Not working)
  for(G4int i=0;i<fNbOfAbsor;i++)
  {
    G4Region* aRegion = new G4Region(layerName[i]);
    fLogicAbsor[i]->SetRegion(aRegion);
    aRegion->AddRootLogicalVolume(fLogicAbsor[i]);
  }
  
  // Attempt Calor Layers as Regions (not working)
  for(G4int i=0;i<fNbOfLayers;i++)
  {
    G4Region* aRegion = new G4Region(layerName[i]);
    fLogicLayer->SetRegion(aRegion);
    aRegion->AddRootLogicalVolume(fLogicLayer);
  }
  */
  // Attempt Whole Calor as One Region
  //G4Region* aRegion = new G4Region(calorRegionName);
  //fLogicCalor->SetRegion(aRegion);
  //aRegion->AddRootLogicalVolume(fLogicCalor);

  // Attempt At Each Layer as the Same Region.
  
    for (G4int j=1; j<=fNbOfAbsor; j++) {
      G4Region* aRegion = new G4Region(absorName[j]);
      fLogicAbsor[j]->SetRegion(aRegion);
      aRegion->AddRootLogicalVolume(fLogicAbsor[j]);
    }
  
  
  
  PrintCalorParameters();

  //always return the fPhysical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The calorimeter is " << fNbOfLayers << " layers of:";
  for (G4int i=1; i<=fNbOfAbsor; i++)
     {
      G4cout << "\n \t" << std::setw(12) << fAbsorMaterial[i]->GetName() <<": "
              << std::setw(6) << G4BestUnit(fAbsorThickness[i],"Length");
     }
  G4cout << "\n-------------------------------------------------------------\n";
  
  G4cout << "\n" << fDefaultMaterial << G4endl;    
  for (G4int j=1; j<=fNbOfAbsor; j++)
     G4cout << "\n" << fAbsorMaterial[j] << G4endl;

  G4cout << "\n-------------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& material)
{
  // search the material by its name
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) fDefaultMaterial = pttoMaterial;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSourceMaterial(const G4String& material)
{
  // search the material by its name
  G4Material* pttoMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) sourceContainerMaterial = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSourceVolumePosition(G4double x , G4double y , G4double z)
{
  // Place the emitting Source at x,y,z... The Placement of SourceVolume will place the emitting source at the proper depth of the volume, with the provided 'windowThickness'
  fSourceVolume_x = x;
  fSourceVolume_y = y;
  fSourceVolume_z = z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCryoMaterial(const G4String& material)
{
  // search the material by its name
  G4Material* pttoMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) cryoMaterial = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCryoCutoutMaterial(const G4String& material)
{
  // search the material by its name
  G4Material* pttoMaterial =
  G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) cryoCutoutMaterial = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfLayers(G4int ival)
{
  // set the number of Layers
  //
  if (ival < 1)
    { G4cout << "\n --->warning from SetfNbOfLayers: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fNbOfLayers = ival;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfAbsor(G4int ival)
{
  // set the number of Absorbers
  //
  if (ival < 1 || ival > (MaxAbsor-1))
    { G4cout << "\n ---> warning from SetfNbOfAbsor: "
             << ival << " must be at least 1 and and most " << MaxAbsor-1
             << ". Command refused" << G4endl;
      return;
    }
  fNbOfAbsor = ival;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorMaterial(G4int ival, const G4String& material)
{
  // search the material by its name
  //
  if (ival > fNbOfAbsor || ival <= 0)
    { G4cout << "\n --->warning from SetAbsorMaterial: absor number "
             << ival << " out of range. Command refused" << G4endl;
      return;
    }

  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) fAbsorMaterial[ival] = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorThickness(G4int ival,G4double val)
{
  // change Absorber thickness
  //
  if (ival > fNbOfAbsor || ival <= 0)
    { G4cout << "\n --->warning from SetAbsorThickness: absor number "
             << ival << " out of range. Command refused" << G4endl;
      return;
    }
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetAbsorThickness: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fAbsorThickness[ival] = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCalorSizeYZ(G4double val)
{
  // change the transverse size
  //
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetfCalorSizeYZ: thickness "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fCalorSizeYZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  //
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(fMagField) delete fMagField;                //delete the existing magn field

  if(fieldValue!=0.)                        // create a new one if non nul
  { fMagField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
    fieldMgr->SetDetectorField(fMagField);
    fieldMgr->CreateChordFinder(fMagField);
  } else {
    fMagField = 0;
    fieldMgr->SetDetectorField(fMagField);
  }
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
