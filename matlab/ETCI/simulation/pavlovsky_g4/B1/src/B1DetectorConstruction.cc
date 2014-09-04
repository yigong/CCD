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
// $Id$
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "B1SteppingAction.hh"
   // use of stepping action to set the accounting volume

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"
#include "G4UnionSolid.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Polycone.hh"

#include "G4PhysicalVolumeStore.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{

  G4double detTemp = 130.0*kelvin;

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // X is depth, Y is height and Z is width
  G4double world_hx  = 0.5*m;
  G4double world_hy  = 0.50*m;
  G4double world_hz  = 0.50*m;

  G4Box* worldBox = new G4Box("World", world_hx, world_hy, world_hz);



  /// ..............................
  /// Start creating geometrical volumes
  /// ..............................



  /// Start table construction
  G4double table_hx = world_hx - 0.1*cm;
  G4double table_hy = 0.5*cm;
  G4double table_hz = world_hz - 0.1*cm;
  G4double table_height_hy = -19.5*cm;

  G4Box* tableBox = new G4Box("Table", table_hx, table_hy, table_hz);

  /// Start source holder construction

  // Definition of source holding block
  G4double sh_on_axis_rot = -M_PI/6*rad;
  G4double sh_height_hx = -0.3*m;
  G4double sh_height_hy = table_height_hy+17.6*cm;
  G4double sh_height_hz = 19.6*cm;
  G4double sh_hole_loc_yz = 4.0*cm;

  G4double sh_hx = 2.7*cm;
  G4double sh_hy = 13.0*cm;
  G4double sh_hz = 13.0*cm;

  G4Box* shBox = new G4Box("Table", sh_hx/2.0, sh_hy/2.0, sh_hz/2.0);

  // Definition of cutouts
  G4double sh_large_cutout_inner_rad = 0.*cm;
  G4double sh_large_cutout_outer_rad = 2.5*cm/2.0;
  G4double sh_large_cutout_hx = sh_hx - 1.0*cm+0.02*cm;

  G4Tubs *shLargeCutoutTubs = new G4Tubs("shLargeCutout",
                                      sh_large_cutout_inner_rad,
                                      sh_large_cutout_outer_rad,
                                      sh_large_cutout_hx/2.0,
                                      0.*deg,
                                      360.*deg);

  G4double sh_small_cutout_inner_rad = 0.*cm;
  G4double sh_small_cutout_outer_rad = 1.2*cm/2.0;
  G4double sh_small_cutout_hx = 1.0*cm+0.02*cm;

  G4Tubs *shSmallCutoutTubs = new G4Tubs("shSmallCutout",
                                      sh_small_cutout_inner_rad,
                                      sh_small_cutout_outer_rad,
                                      sh_small_cutout_hx/2.0,
                                      0.*deg,
                                      360.*deg);

  G4ThreeVector xTransSH_1(0.5*(-sh_large_cutout_hx+sh_hx)+0.01*cm, sh_hole_loc_yz, sh_hole_loc_yz);
  G4RotationMatrix* yRotSH = new G4RotationMatrix();
  yRotSH->rotateY(M_PI/2.0*rad);

  G4SubtractionSolid* shBoxSub_1 = new G4SubtractionSolid("shBoxSub_1", shBox, shLargeCutoutTubs, yRotSH, xTransSH_1);

  G4ThreeVector xTransSH(-0.5*(-sh_small_cutout_hx+sh_hx)-0.01*cm, sh_hole_loc_yz, sh_hole_loc_yz);
  G4SubtractionSolid* shBoxSub_2 = new G4SubtractionSolid("shBoxSub_2", shBoxSub_1, shSmallCutoutTubs, yRotSH, xTransSH);



  xTransSH_1.set(0.5*(-sh_large_cutout_hx+sh_hx)+0.01*cm, 1.0/3.0*sh_hole_loc_yz, 1.0/3.*sh_hole_loc_yz);
  G4SubtractionSolid* shBoxSub_3 = new G4SubtractionSolid("shBoxSub_3", shBoxSub_2, shLargeCutoutTubs, yRotSH, xTransSH_1);

  xTransSH.set(-0.5*(-sh_small_cutout_hx+sh_hx)-0.01*cm, 1./3.*sh_hole_loc_yz, 1./3.*sh_hole_loc_yz);
  G4SubtractionSolid* shBoxSub_4 = new G4SubtractionSolid("shBoxSub", shBoxSub_3, shSmallCutoutTubs, yRotSH, xTransSH);


  xTransSH_1.set(0.5*(-sh_large_cutout_hx+sh_hx)+0.01*cm, -1.0/3.0*sh_hole_loc_yz, -1.0/3.*sh_hole_loc_yz);
  G4SubtractionSolid* shBoxSub_5 = new G4SubtractionSolid("shBoxSub_3", shBoxSub_4, shLargeCutoutTubs, yRotSH, xTransSH_1);

  xTransSH.set(-0.5*(-sh_small_cutout_hx+sh_hx)-0.01*cm, -1./3.*sh_hole_loc_yz, -1./3.*sh_hole_loc_yz);
  G4SubtractionSolid* shBoxSub_6 = new G4SubtractionSolid("shBoxSub", shBoxSub_5, shSmallCutoutTubs, yRotSH, xTransSH);


  xTransSH_1.set(0.5*(-sh_large_cutout_hx+sh_hx)+0.01*cm, -sh_hole_loc_yz, -sh_hole_loc_yz);
  G4SubtractionSolid* shBoxSub_7 = new G4SubtractionSolid("shBoxSub_3", shBoxSub_6, shLargeCutoutTubs, yRotSH, xTransSH_1);

  xTransSH.set(-0.5*(-sh_small_cutout_hx+sh_hx)-0.01*cm, -sh_hole_loc_yz, -sh_hole_loc_yz);
  G4SubtractionSolid* shBoxSub = new G4SubtractionSolid("shBoxSub", shBoxSub_7, shSmallCutoutTubs, yRotSH, xTransSH);

  //25 12 40


  /// Start Cryostat definitions

  // Definition of cryostat external dims
  G4double cryo_hx  = 6.1*cm;
  G4double cryo_hy  = 17.3*cm;
  G4double cryo_hz  = 23.4*cm;

  G4Box* cryoBox = new G4Box("Cryo", cryo_hx/2.0, cryo_hy/2.0, cryo_hz/2.0);

  // Definition of vacuum volume
  G4double cryo_thickness_x = 0.5*cm;
  G4double cryo_thickness_y = 0.7*cm;
  G4double cryo_thickness_z = cryo_thickness_y;

  G4double cryo_cutout_hx = cryo_hx-2*cryo_thickness_x;
  G4double cryo_cutout_hy = cryo_hy-2*cryo_thickness_y;
  G4double cryo_cutout_hz = cryo_hz-2*cryo_thickness_z;

  G4Box* cryoCutoutBox = new G4Box("CryoCutout", cryo_cutout_hx/2.0, cryo_cutout_hy/2.0, cryo_cutout_hz/2.0);

  // Definition of cryostat peripherial attachments
  G4double top_plate_hx = 5.50*cm;
  G4double top_plate_hy = 0.6*cm;
  G4double top_plate_hz = 16.26*cm;

  G4Box* topPlateBox = new G4Box("TopPlate", top_plate_hx/2.0, top_plate_hy/2.0, top_plate_hz/2.0);

  // Definition of front cutout
  G4double front_cutout_inner_rad = 0.*cm;
  G4double front_cutout_outer_rad = 4.99*cm;
  G4double front_cutout_hx = 0.4*cm;

  G4Tubs *frontCutoutTubs = new G4Tubs("frontCutout",
                                      front_cutout_inner_rad,
                                      front_cutout_outer_rad,
                                      front_cutout_hx/2.0,
                                      0.*deg,
                                      360.*deg);

  G4ThreeVector xTransFrontCutout_1(-0.5*(cryo_hx-front_cutout_hx)-0.01*cm, 0.*cm, 0.*cm);
  G4RotationMatrix* yRotFrontCutout = new G4RotationMatrix();
  yRotFrontCutout->rotateY(M_PI/2.0*rad);

  G4SubtractionSolid* cryoBoxSub_1 = new G4SubtractionSolid("CyroCutoutComplete", cryoBox, frontCutoutTubs, yRotFrontCutout, xTransFrontCutout_1);

  G4ThreeVector xTransFrontCutout(0.5*(cryo_hx-front_cutout_hx)+0.01*cm, 0.*cm, 0.*cm);
  G4SubtractionSolid* cryoBoxSub = new G4SubtractionSolid("CyroCutoutComplete", cryoBoxSub_1, frontCutoutTubs, yRotFrontCutout, xTransFrontCutout);

  // Definition of vacuum port(s)
  G4double vac_port_inner_rad = 0.6*cm;
  G4double vac_port_outer_rad = 1.921/2.0*cm;
  G4double vac_port_hy = 2.0*cm;

  G4Tubs *vacPortTubs = new G4Tubs("vacPort",
                                      vac_port_inner_rad,
                                      vac_port_outer_rad,
                                      vac_port_hy/2.0,
                                      0.*deg,
                                      360.*deg);

  G4double vac_port_inner_rad_2 = 0.6*cm;
  G4double vac_port_outer_rad_2 = 3.0/2.0*cm;
  G4double vac_port_hy_2 = 2.4*cm;

  G4Tubs *vacPortTubs_2 = new G4Tubs("vacPort_2",
                                      vac_port_inner_rad_2,
                                      vac_port_outer_rad_2,
                                      vac_port_hy_2/2.0,
                                      0.*deg,
                                      360.*deg);

  G4double vac_port_inner_rad_3 = 3.0/2.0*cm;
  G4double vac_port_outer_rad_3 = 3.3/2.0*cm;
  G4double vac_port_hy_3 = 5.3*cm;

  G4Tubs *vacPortTubs_3 = new G4Tubs("vacPort_3",
                                      vac_port_inner_rad_3,
                                      vac_port_outer_rad_3,
                                      vac_port_hy_3/2.0,
                                      0.*deg,
                                      360.*deg);

  /// Start electronics boxes defintion

  // Define connector box
  G4double connector_hx  = 5.1*cm;
  G4double connector_hy  = 16.3*cm;
  G4double connector_hz  = 2.526*cm;

  G4Box* connectorBox = new G4Box("ConnectorBox", connector_hx/2.0, connector_hy/2.0, connector_hz/2.0);

  // Define connector box cutout
  G4double connector_thickness = 0.28*cm;
  G4double connector_cutout_hx  = connector_hx-2.0*connector_thickness;
  G4double connector_cutout_hy  = connector_hy-2.0*connector_thickness;
  G4double connector_cutout_hz  = connector_hz-2.0*connector_thickness;

  G4Box* connectorCutoutBox = new G4Box("ConnectorCutoutBox", connector_cutout_hx/2.0, connector_cutout_hy/2.0, connector_cutout_hz/2.0);

  // Define interface board box
  G4double if_hx  = 8.3*cm;
  G4double if_hy  = 19.5*cm;
  G4double if_hz  = 6.0*cm;

  G4Box* ifBox = new G4Box("IFBox", if_hx/2.0, if_hy/2.0, if_hz/2.0);

  // Define interface board box
  G4double if_thickness = 0.28*cm;
  G4double if_cutout_hx  = if_hx - 2.0*if_thickness;
  G4double if_cutout_hy  = if_hy - 2.0*if_thickness;
  G4double if_cutout_hz  = if_hz - 2.0*if_thickness;

  G4Box* ifCutoutBox = new G4Box("IFCutoutBox", if_cutout_hx/2.0, if_cutout_hy/2.0, if_cutout_hz/2.0);


  /// Start dewar construction

  // Definition of dewar port
  G4double dewar_port_inner_rad = 0;
  G4double dewar_port_outer_rad = 4.421/2.0*cm;
  G4double dewar_port_hy = 7.851*cm;

  G4Tubs *dewarPortTubs = new G4Tubs("dewarPort",
                                      dewar_port_inner_rad,
                                      dewar_port_outer_rad,
                                      dewar_port_hy/2.0,
                                      0.*deg,
                                      360.*deg);

  // Definition of dewar port cutout
  G4double dewar_port_cutout_inner_rad = 0;
  G4double dewar_port_cutout_outer_rad = dewar_port_outer_rad-0.2*cm;
  G4double dewar_port_cutout_hy = dewar_port_hy-0.2*cm;

  G4Tubs *dewarPortCutoutTubs = new G4Tubs("dewarPortCutout",
                                      dewar_port_cutout_inner_rad,
                                      dewar_port_cutout_outer_rad,
                                      dewar_port_cutout_hy/2.0,
                                      0.*deg,
                                      360.*deg);

  // Definition of dewar
  G4double dewar_inner_rad = 0;
  G4double dewar_outer_rad = 21.5/2.0*cm;
  G4double dewar_hy = 17.5*2.54*cm;
  G4double z_locs[] = {-dewar_hy/2.0, dewar_hy/2.0-5.5*cm, dewar_hy/2.0};
  G4double r_inner_values[] = {0., 0.,0.};
  G4double r_outer_values[] = {dewar_outer_rad, dewar_outer_rad, 4.*cm};

  G4Polycone* dewarPoly = new G4Polycone("dewar",0., M_PI*2.0*rad, 3., z_locs, r_inner_values, r_outer_values);


  // Definition of dewar cutout
  G4double dewar_cutout_inner_rad = 0;
  G4double dewar_cutout_outer_rad = dewar_outer_rad-0.2*cm;
  G4double dewar_cutout_hy = dewar_hy-0.2*cm;
  G4double z_cutout_locs[] = {-dewar_cutout_hy/2.0, dewar_cutout_hy/2.0-5.5*cm, dewar_cutout_hy/2.0};
  G4double r_cutout_inner_values[] = {0., 0.,0.};
  G4double r_cutout_outer_values[] = {dewar_cutout_outer_rad, dewar_cutout_outer_rad, 4.*cm-0.2*cm};

  G4Polycone* dewarCutoutPoly = new G4Polycone("dewarCutout",0., M_PI*2.0*rad, 3., z_cutout_locs, r_cutout_inner_values, r_cutout_outer_values);

  // Definition of liquid nitrogen volume
  G4double ln_vol_inner_rad = 0;
  G4double ln_vol_outer_rad = dewar_outer_rad-2.*cm;
  G4double ln_vol_hy = dewar_hy-9*cm;

  G4Tubs *lnVolTubs = new G4Tubs("lnVol",
                                      ln_vol_inner_rad,
                                      ln_vol_outer_rad,
                                      ln_vol_hy/2.0,
                                      0.*deg,
                                      360.*deg);


  /// Start CCD definitions

  // CCD definition
  G4double ccd_hx  = 650.*um;
  G4double ccd_hy  = 1.5*cm;
  G4double ccd_hz  = 0.75*cm;

  G4Box* ccdBox = new G4Box("ccd", ccd_hx/2.0, ccd_hy/2.0, ccd_hz/2.0);


  // Definition of cryostat external dims
  G4double ccd_mount_board_hx  = 0.15*cm;
  G4double ccd_mount_board_hy  = 8.*cm;
  G4double ccd_mount_board_hz  = 7.5*cm;

  G4Box* ccdMountBoardBox = new G4Box("ccdMountBoard", ccd_mount_board_hx/2.0, ccd_mount_board_hy/2.0, ccd_mount_board_hz/2.0);


  // CCD definition
  G4double ccd_cutout_hx  = 2.5*650*um;
  G4double ccd_cutout_hy  = 1.1*1.5*cm;
  G4double ccd_cutout_hz  = 1.2*0.75*cm;

  G4Box* ccdCutoutBox = new G4Box("ccd", ccd_cutout_hx/2.0, ccd_cutout_hy/2.0, ccd_cutout_hz/2.0);

  // Mounting board definition
  G4SubtractionSolid* ccdMountBoardBoxSub = new G4SubtractionSolid("ccdMountingBoardBoxSub", ccdMountBoardBox, ccdCutoutBox);

  // Define aluminum nitride backing
  G4double aln_hx  = 0.5*mm;
  G4double aln_hy  = ccd_mount_board_hy;
  G4double aln_hz  = 1./3.*ccd_mount_board_hz;

  G4Box* alNitrideBackBox = new G4Box("AlNBack", aln_hx/2.0, aln_hy/2.0, aln_hz/2.0);

  // Define picture frame board
  G4double pfb_hx  = 1.*mm;
  G4double pfb_hy  = 76.21*mm;
  G4double pfb_hz  = 101.67*mm;
  G4double pfb_cold_mass_offset_hx = 1.*cm;
  G4double pfb_cold_mass_offset_hy = 1.5*cm;

  G4Box* pfbBox = new G4Box("pfb", pfb_hx/2.0, pfb_hy/2.0, pfb_hz/2.0);

  // Define picture frame board cutout
  G4double pfb_cutout_hx  = 1.1*mm;
  G4double pfb_cutout_hy  = 39.52*mm;
  G4double pfb_cutout_hz  = 38.20*mm;
  G4double pfb_cutout_offset_y = 1.*cm;

  G4ThreeVector xTransPFB(0.*cm, (pfb_cutout_offset_y - (pfb_hy-pfb_cutout_hy)/2.), 0.*cm);
  G4RotationMatrix* yRotPFB = new G4RotationMatrix;

  G4Box* pfbCutoutBox = new G4Box("pictureFrameBoard", pfb_cutout_hx/2.0, pfb_cutout_hy/2.0, pfb_cutout_hz/2.0);

  G4SubtractionSolid* pfbSubBox = new G4SubtractionSolid("pfb", pfbBox, pfbCutoutBox,yRotPFB,xTransPFB);


  /// Start Cold plate definitions

  // coldplate
  G4double cold_hx  = 0.5*cm;
  G4double cold_hy  = 12.425*cm;
  G4double cold_hz  = 9.9930*cm;

  G4Box* coldBox = new G4Box("coldplate", cold_hx/2.0, cold_hy/2.0, cold_hz/2.0);

  // cold plate cut
  G4double cold_cutout_hx  = cold_hx+0.1*cm;
  G4double cold_cutout_hy  = 6.4178*cm;
  G4double cold_cutout_hz  = 7.141*cm;
  G4double cold_cutout_offset_hy = 2.122*cm;

  G4Box* coldCutoutBox = new G4Box("coldplatecutout", cold_cutout_hx/2.0, cold_cutout_hy/2.0, cold_cutout_hz/2.0);

  // make cold finger clamp
  G4double cold_clamp_hx = 2.0*cm;
  G4double cold_clamp_hy = 5.0*cm;
  G4double cold_clamp_hz = 3.0*cm;

  G4Box* coldClampBox = new G4Box("coldclamp", cold_clamp_hx/2.0, cold_clamp_hy/2.0, cold_clamp_hz/2.0);

  // make cold finger
  G4double cold_finger_inner_rad = 0.*cm;
  G4double cold_finger_outer_rad = 0.75*cm;
  G4double cold_finger_hz = 4.*cm;
  G4double cold_finger_offset = 1.*cm;

  G4Tubs *coldFingerTubs = new G4Tubs("coldFinger",
                                      cold_finger_inner_rad,
                                      cold_finger_outer_rad,
                                      cold_finger_hz/2.0,
                                      0.*deg,
                                      360.*deg);

  // make cold finger clamp with cutout
  G4ThreeVector xTransColdFinger(0.*cm, 0.*cm, 0.5*(cold_clamp_hz-cold_finger_hz) - cold_finger_offset);
  G4RotationMatrix* yRotColdFinger = new G4RotationMatrix;

  G4SubtractionSolid* coldClampSub = new G4SubtractionSolid("coldClamp", coldClampBox, coldFingerTubs, yRotColdFinger, xTransColdFinger );

  // make cold plate
  G4ThreeVector xTransVec(0.*cm, 0.5*(cold_hy-cold_cutout_hy) - cold_cutout_offset_hy, 0.5*(cold_hz-cold_cutout_hz));
  G4RotationMatrix* yRot = new G4RotationMatrix;

  G4SubtractionSolid* coldPlateSub = new G4SubtractionSolid("coldPlate", coldBox, coldCutoutBox, yRot, xTransVec);

  // make copper bar
  G4double copper_bar_hx  = 1.5*cm;
  G4double copper_bar_hy  = 1.*cm;
  G4double copper_bar_hz  = 9.6*cm;

  G4Box* copperBarBox = new G4Box("copperBar", copper_bar_hx/2.0, copper_bar_hy/2.0, copper_bar_hz/2.0);



  /// ..............................
  /// Start loading materials
  /// ..............................



  G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* Al = nist->FindOrBuildMaterial("G4_Al");
//  G4Material* Al = nist->FindOrBuildMaterial("G4_Pb");
  G4Material* Vac = new G4Material("Vac", 14., 28.*g/mole, 1.e-7*mole/m3,
                     kStateGas, detTemp, 1.e-6*pascal);
  G4Material* Si = nist->FindOrBuildMaterial("G4_Si");
  G4Material* Cu = nist->FindOrBuildMaterial("G4_Cu");
  G4Material* LN = nist->FindOrBuildMaterial("G4_lN2");

  G4Element* eAl = nist->FindOrBuildElement("Al");
  G4Element* eN = nist->FindOrBuildElement("N");

  G4Material* AlN = new G4Material("AlN", 3.260*g/cm3 ,2);
  AlN->AddElement(eN, 1);
  AlN->AddElement(eAl, 1);



  /// ..............................
  /// Start creating logical volumes
  /// ..............................



  G4LogicalVolume* worldLog = new G4LogicalVolume(worldBox, Air, "World");
  G4LogicalVolume* tableLog = new G4LogicalVolume(tableBox, Al, "Table");
  G4LogicalVolume* shLog = new G4LogicalVolume(shBoxSub, Al, "SourceHolder");
  G4LogicalVolume* cryoLog = new G4LogicalVolume(cryoBoxSub, Al, "Cryo");
  G4LogicalVolume* vacPortLog = new G4LogicalVolume(vacPortTubs, Al, "VacPort");
  G4LogicalVolume* vacPortLog_2 = new G4LogicalVolume(vacPortTubs_2, Al, "VacPort_2");
  G4LogicalVolume* vacPortLog_3 = new G4LogicalVolume(vacPortTubs_3, Al, "VacPort_3");
  G4LogicalVolume* cryoCutoutLog = new G4LogicalVolume(cryoCutoutBox, Vac, "CryoCutout");
  G4LogicalVolume* topPlateLog = new G4LogicalVolume(topPlateBox, Al, "TopPlate");
  G4LogicalVolume* connectorLog = new G4LogicalVolume(connectorBox, Al, "ConnectorBox");
  G4LogicalVolume* connectorCutoutLog = new G4LogicalVolume(connectorCutoutBox, Air, "ConnectorBoxInternal");
  G4LogicalVolume* ifLog = new G4LogicalVolume(ifBox, Al, "IFBox");
  G4LogicalVolume* ifCutoutLog = new G4LogicalVolume(ifCutoutBox, Air, "IFBoxInternal");
  G4LogicalVolume* dewarPortLog = new G4LogicalVolume(dewarPortTubs, Al, "DewarPort");
  G4LogicalVolume* dewarPortCutoutLog = new G4LogicalVolume(dewarPortCutoutTubs, Vac, "DewarPortCutout");
  G4LogicalVolume* dewarLog = new G4LogicalVolume(dewarPoly, Al, "Dewar");
  G4LogicalVolume* dewarCutoutLog = new G4LogicalVolume(dewarCutoutPoly, Vac, "DewarCutout");
  G4LogicalVolume* lnVolLog = new G4LogicalVolume(lnVolTubs, LN, "LnVol");
  G4LogicalVolume* ccdLog = new G4LogicalVolume(ccdBox, Si, "Ccd");
  G4LogicalVolume* ccdMountBoardLog = new G4LogicalVolume(ccdMountBoardBoxSub, Si, "MountingBoard");
  G4LogicalVolume* alNitrideBackLog = new G4LogicalVolume(alNitrideBackBox, AlN, "AlNBacking");
  G4LogicalVolume* pfbLog = new G4LogicalVolume(pfbSubBox, Si, "pfb");
  G4LogicalVolume* coldClampLog = new G4LogicalVolume(coldClampSub, Al, "ColdClamp");
  G4LogicalVolume* coldPlateLog = new G4LogicalVolume(coldPlateSub, Al, "ColdPlate");
  G4LogicalVolume* coldFingerLog = new G4LogicalVolume(coldFingerTubs, Cu, "ColdPlate");
  G4LogicalVolume* copperBarLog = new G4LogicalVolume(copperBarBox, Cu, "CopperBar");

  G4Region *ccdLogReg = new G4Region("Ccd");
  ccdLogReg->AddRootLogicalVolume(ccdLog);

  G4Region *ccdMountBoardLogReg = new G4Region("ccdMountBoard");
  ccdMountBoardLogReg->AddRootLogicalVolume(ccdMountBoardLog);

  G4Region *alNitrideBackLogReg = new G4Region("alNitrideBack");
  alNitrideBackLogReg->AddRootLogicalVolume(alNitrideBackLog);

  G4Region *pfbLogReg = new G4Region("pfb");
  pfbLogReg->AddRootLogicalVolume(pfbLog);

  G4Region *coldPlateLogReg = new G4Region("coldPlate");
  coldPlateLogReg->AddRootLogicalVolume(coldPlateLog);

  G4Region *copperBarLogReg = new G4Region("copperBar");
  copperBarLogReg->AddRootLogicalVolume(copperBarLog);

  G4Region *cryoLogReg = new G4Region("cryoLog");
  cryoLogReg->AddRootLogicalVolume(cryoLog);

  G4VisAttributes * lnVolColor = new G4VisAttributes(G4Colour(0.,1.,1.));
  // Set the forced wireframe style
  lnVolColor->SetForceSolid(true);
  // Assignment of the visualization attributes to the logical volume
  lnVolLog->SetVisAttributes(lnVolColor);


  G4VisAttributes * copperColor = new G4VisAttributes(G4Colour::Yellow());
  // Set the forced wireframe style
  copperColor->SetForceSolid(true);
  // Assignment of the visualization attributes to the logical volume
  copperBarLog->SetVisAttributes(copperColor);
  coldFingerLog->SetVisAttributes(copperColor);

  G4VisAttributes * alColor = new G4VisAttributes(G4Colour::Grey());
  // Set the forced wireframe style
  alColor->SetForceSolid(true);
  // Assignment of the visualization attributes to the logical volume
  coldPlateLog->SetVisAttributes(alColor);
  coldClampLog->SetVisAttributes(alColor);
  shLog->SetVisAttributes(alColor);
  tableLog->SetVisAttributes(alColor);

  G4VisAttributes * ccdColor = new G4VisAttributes(G4Colour::Red());
  ccdColor->SetForceSolid(true);
  ccdLog->SetVisAttributes(ccdColor);

  G4VisAttributes * ccdMountColor = new G4VisAttributes(G4Colour::Blue());
  ccdMountColor->SetForceSolid(true);
  ccdMountBoardLog->SetVisAttributes(ccdMountColor);

  G4VisAttributes * alnColor = new G4VisAttributes(G4Colour::Magenta());
  alnColor->SetForceSolid(true);
  alNitrideBackLog->SetVisAttributes(alnColor);

  G4VisAttributes * siColor = new G4VisAttributes(G4Colour::Green());
  siColor->SetForceSolid(true);
  pfbLog->SetVisAttributes(siColor);


  /// ..............................
  /// Start creating phsyical volumes
  /// ..............................



  G4double pos_x = 0.*world_hx;
  G4double pos_y = 0.*world_hy;
  G4double pos_z = 0.*world_hz;

  G4double end_space_z = 0.5*cm;

  G4double cold_mass_offset_z = 0.07*cold_hz;
  G4double cold_mass_offset_y = 0.07*cold_hy;

  G4VPhysicalVolume* cryoPhys = new G4PVPlacement(0, //no rotation
                                                     G4ThreeVector(pos_x,pos_y,pos_z),
                                                     cryoLog,
                                                     "Cryo",
                                                     worldLog,
                                                     false,
                                                     0);

  G4VPhysicalVolume* cryoCutoutPhys = new G4PVPlacement(0,
                                                        G4ThreeVector(0,0,0),
                                                        cryoCutoutLog,
                                                        "CryoCutout",
                                                        cryoLog,
                                                        false,
                                                        0);

//  G4VPhysicalVolume* topPlatePhys = new G4PVPlacement(0,
//                                                        G4ThreeVector(pos_x,
//                                                                      pos_y+0.5*(cryo_hy+top_plate_hy),
//                                                                      pos_z-end_space_z+0.5*(cryo_hz-top_plate_hz)),
//                                                        topPlateLog,
//                                                        "TopPlate",
//                                                        worldLog,
//                                                        false,
//                                                        0);

//  G4RotationMatrix* rotVacPort = new G4RotationMatrix();
//  rotVacPort->rotateX(M_PI/2.0*rad);

//  G4VPhysicalVolume* vacPortPhys = new G4PVPlacement(rotVacPort,
//                                                        G4ThreeVector(pos_x,
//                                                                      pos_y+0.5*(cryo_hy+vac_port_hy),
//                                                                      pos_z-0.5*cryo_hz+2.5*cm),
//                                                        vacPortLog,
//                                                        "VacPort",
//                                                        worldLog,
//                                                        false,
//                                                        0);

//  G4VPhysicalVolume* vacPortPhys_2 = new G4PVPlacement(rotVacPort,
//                                                        G4ThreeVector(pos_x,
//                                                                      pos_y+0.5*(cryo_hy+vac_port_hy_2)+vac_port_hy,
//                                                                      pos_z-0.5*cryo_hz+2.5*cm),
//                                                        vacPortLog_2,
//                                                        "VacPort_2",
//                                                        worldLog,
//                                                        false,
//                                                        0);

//  G4VPhysicalVolume* vacPortPhys_3 = new G4PVPlacement(rotVacPort,
//                                                        G4ThreeVector(pos_x,
//                                                                      pos_y-0.5*(cryo_hy+vac_port_hy_3),
//                                                                      pos_z+0.5*cryo_hz-2.5*cm),
//                                                        vacPortLog_3,
//                                                        "VacPort_3",
//                                                        worldLog,
//                                                        false,
//                                                        0);

//  G4VPhysicalVolume* connectorBoxPhys = new G4PVPlacement(0, //no rotation
//                                                  G4ThreeVector(pos_x,pos_y,pos_z+0.5*(cryo_hz+connector_hz)),
//                                                     connectorLog,
//                                                     "ConnectorBox",
//                                                     worldLog,
//                                                     false,
//                                                     0);


//  G4VPhysicalVolume* connectorCutoutPhys = new G4PVPlacement(0, //no rotation
//                                                  G4ThreeVector(0.,0.,0.),
//                                                     connectorCutoutLog,
//                                                     "ConnectorCutout",
//                                                     connectorLog,
//                                                     false,
//                                                     0);

//  G4VPhysicalVolume* ifBoxPhys = new G4PVPlacement(0, //no rotation
//                                                   G4ThreeVector(pos_x,pos_y-0.5*(if_hy-connector_hy)+0.28*cm,pos_z+0.5*(cryo_hz+if_hz)+connector_hz),
//                                                     ifLog,
//                                                     "ConnectorBox",
//                                                     worldLog,
//                                                     false,
//                                                     0);


//  G4VPhysicalVolume* ifCutoutPhys = new G4PVPlacement(0, //no rotation
//                                                  G4ThreeVector(0.,0.,0.),
//                                                     ifCutoutLog,
//                                                     "ConnectorCutout",
//                                                     ifLog,
//                                                     false,
//                                                     0);

//  G4VPhysicalVolume* dewarPortPhys = new G4PVPlacement(0,
//                                                  G4ThreeVector(0,-cold_mass_offset_y, -0.5*(cryo_hz+dewar_port_hy)),
//                                                     dewarPortLog,
//                                                     "DewarPort",
//                                                     worldLog,
//                                                     false,
//                                                     0);

//  G4VPhysicalVolume* dewarPortCutoutPhys = new G4PVPlacement(0, //no rotation
//                                                  G4ThreeVector(0.,0.,0.),
//                                                     dewarPortCutoutLog,
//                                                     "DewarPortCutout",
//                                                     dewarPortLog,
//                                                     false,
//                                                     0);

//  G4VPhysicalVolume* dewarPhys = new G4PVPlacement(rotVacPort, //no rotation
//                                                   G4ThreeVector(0.,0.5*(dewar_hy-cryo_hy)+1.*cm,-0.5*(cryo_hz+2.*dewar_outer_rad)-dewar_port_hy),
//                                                     dewarLog,
//                                                     "dewar",
//                                                     worldLog,
//                                                     false,
//                                                     0);

//  G4VPhysicalVolume* dewarCutoutPhys = new G4PVPlacement(0, //no rotation
//                                                  G4ThreeVector(0.,0.,0.),
//                                                     dewarCutoutLog,
//                                                     "dewarCutout",
//                                                     dewarLog,
//                                                     false,
//                                                     0);

//  G4VPhysicalVolume* lnVolPhys = new G4PVPlacement(0,
//                                                   G4ThreeVector(0.,0.,2.5*cm-0.5*(dewar_cutout_hy-ln_vol_hy)),
//                                                   lnVolLog,
//                                                   "lnvol",
//                                                   dewarCutoutLog,
//                                                   false,
//                                                   0);

  G4VPhysicalVolume* ccdMountBoardPhys = new G4PVPlacement(0,
                                                           G4ThreeVector(0.5*(ccd_mount_board_hx+cold_hx),0,0),
                                                 ccdMountBoardLog,
                                                 "CcdMountingBoard",
                                                 cryoCutoutLog,
                                                 false,
                                                 0);

  G4VPhysicalVolume* alNitrideBackPhys = new G4PVPlacement(0,
                                                 G4ThreeVector(0.5*(aln_hx+cold_hx)+ccd_mount_board_hx,0,0),
                                                 alNitrideBackLog,
                                                 "CcdMountingBoard",
                                                 cryoCutoutLog,
                                                 false,
                                                 0);


  G4VPhysicalVolume* ccdPhys = new G4PVPlacement(0,
                                                 G4ThreeVector(0.5*(cold_hx-ccd_hx)+ccd_mount_board_hx,0,0),
                                                 ccdLog,
                                                 "Ccd",
                                                 cryoCutoutLog,
                                                 false,
                                                 0);

  G4VPhysicalVolume* pfbPhys = new G4PVPlacement(0,
                                                 G4ThreeVector(0.5*(aln_hx+cold_hx+ccd_mount_board_hx)+pfb_hx+pfb_cold_mass_offset_hx,
                                                               -(pfb_cold_mass_offset_hy-0.5*(cold_hy-pfb_hy)),0),
                                                 pfbLog,
                                                 "PFB",
                                                 cryoCutoutLog,
                                                 false,
                                                 0);

  G4VPhysicalVolume* coldPlatePhys = new G4PVPlacement(0,
                                                 G4ThreeVector(0,-cold_mass_offset_y,-cold_mass_offset_z),
                                                 coldPlateLog,
                                                 "ColdPlate",
                                                 cryoCutoutLog,
                                                 false,
                                                 0);

  G4VPhysicalVolume* coldFingerPhys = new G4PVPlacement(0,
                                                 G4ThreeVector(0,-cold_mass_offset_y,0.5*(-cold_finger_hz-cold_hz)-cold_mass_offset_z-cold_finger_offset),
                                                 coldFingerLog,
                                                 "ColdFinger",
                                                 cryoCutoutLog,
                                                 false,
                                                 0);

  G4VPhysicalVolume* copperBarPhys = new G4PVPlacement(0,
                                                 G4ThreeVector(0.5*(cold_hx+copper_bar_hx)+ccd_mount_board_hx+aln_hx,
                                                               -0.5*(cold_hy+copper_bar_hy) + cold_cutout_offset_hy + cold_mass_offset_y,
                                                               -cold_mass_offset_z),
                                                 copperBarLog,
                                                 "ColdFinger",
                                                 cryoCutoutLog,
                                                 false,
                                                 0);


  G4VPhysicalVolume* coldClampPhys = new G4PVPlacement(0,
                                                 G4ThreeVector(0,-cold_mass_offset_y,0.5*(-cold_clamp_hz-cold_hz)-cold_mass_offset_z),
                                                 coldClampLog,
                                                 "ColdFinger",
                                                 cryoCutoutLog,
                                                 false,
                                                 0);


//  G4VPhysicalVolume* tablePhys = new G4PVPlacement(0,
//                                                    G4ThreeVector(0,table_height_hy,0),
//                                                    tableLog,
//                                                    "World",
//                                                    worldLog,
//                                                    false,
//                                                    0);

//  G4RotationMatrix* rotSH = new G4RotationMatrix();
//  rotSH->rotateY(sh_on_axis_rot);
//  G4VPhysicalVolume* shPhys = new G4PVPlacement(rotSH,
//                                                    G4ThreeVector(sh_height_hx,sh_height_hy,sh_height_hz),
//                                                    shLog,
//                                                    "SourceHolder",
//                                                    worldLog,
//                                                    false,
//                                                    0);

  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,
                                                    G4ThreeVector(0,0,0),
                                                    worldLog,
                                                    "World",
                                                    NULL,
                                                    false,
                                                    0);


  B1SteppingAction* steppingAction = B1SteppingAction::Instance();
  ////steppingAction->SetVolume(logicShape1);
  steppingAction->SetVolume(ccdLog);


  //
  //always return the physical World
  //

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

