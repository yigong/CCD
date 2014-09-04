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
/// \file electromagnetic/TestEm3/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4Tubs;
class G4SubtractionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;

     const G4int MaxAbsor = 10;                        // 0 + 9  
     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
   DetectorConstruction();
  ~DetectorConstruction();

public:
  
  void SetNbOfAbsor     (G4int);      
  void SetAbsorMaterial (G4int,const G4String&);     
  void SetAbsorThickness(G4int,G4double);
          
  void SetWorldMaterial (const G4String&);
  void SetCalorSizeYZ   (G4double);          
  void SetNbOfLayers    (G4int);
  
  // Set Source Voulme
  void SetSourceMaterial(const G4String&);
  void SetSourceVolumePosition(G4double, G4double, G4double);

  // Set Materials
  void SetCryoMaterial(const G4String&);
  void SetCryoCutoutMaterial(const G4String&);
    
  //void SetGeSizeYZ(G4double);
    
  //void SetMagField   (G4double);
  
  virtual   
  G4VPhysicalVolume* Construct();

  void UpdateGeometry();
     
public:
  
  void PrintCalorParameters(); 
                    
  G4double GetWorldSizeX()            {return fWorldSizeX;};
  G4double GetWorldSizeYZ()           {return fWorldSizeYZ;};
     
  G4double GetCalorThickness()        {return fCalorThickness;};
  G4double GetCalorSizeYZ()           {return fCalorSizeYZ;};
  G4double GetSourceWindowThickness() {return fSourceWindowThickness;};
  
  // New for Cryo and Catcher Detector -------------
    //  G4double GetGeSizeYZ()		{return GeSizeYZ;};
    //  G4double    GetGeThickness()       {return GeThickness;};
    G4double GetCalorOffset()           {return CalorOffset;};
    G4Material* SetCryoCutoutMaterial()	{return cryoCutoutMaterial;};
    G4Material* GetCryoMaterial()       {return cryoMaterial;};
  
    
  //--------------
  G4int GetNbOfLayers()              {return fNbOfLayers;};
     
  G4int       GetNbOfAbsor()             {return fNbOfAbsor;}; 
  G4Material* GetAbsorMaterial(G4int i)  {return fAbsorMaterial[i];};
  G4double    GetAbsorThickness(G4int i) {return fAbsorThickness[i];};      
  // World
  const G4VPhysicalVolume* GetphysiWorld()        {return fPhysiWorld;};
  const G4Material*        GetWorldMaterial()     {return fDefaultMaterial;};
  // Source
  const G4VPhysicalVolume* GetPhysiSource()       {return fPhysiSource;};
  const G4Material*        GetSourceMaterial()    {return sourceContainerMaterial;};

  // Absorber
  const G4VPhysicalVolume* GetAbsorber(G4int i)   {return fPhysiAbsor[i];};
  // Cryo
  // Calor
  
  //----------------

private:
  // Calrimeter Parameters
  G4int              fNbOfAbsor;
  //G4String           layerName[MaxAbsor];
  G4String           calorRegionName;
  G4Material*        fAbsorMaterial [MaxAbsor];
  G4double           fAbsorThickness[MaxAbsor];

  G4int              fNbOfLayers;
  G4double           fLayerThickness;

  G4double           fCalorSizeYZ;
  G4double           fCalorThickness;
  G4double           fCalorOffset;
  
  // X is depth, Y is height and Z is width
  // Can be used to Shift the center of the Cryostat off of 0,0,0
  G4double           world_hx;
  G4double           world_hy;
  G4double           world_hz;
  
  // Definition of cryostat external dims
  G4double           cryo_hx;
  G4double           cryo_hy;
  G4double           cryo_hz;

  // Definition of vacuum volume
  G4double           cryo_thickness_x;
  G4double           cryo_thickness_y;
  G4double           cryo_thickness_z;
  
  // Definition of vacuum volume
  G4double           cryo_cutout_hx;
  G4double           cryo_cutout_hy;
  G4double           cryo_cutout_hz;
  
  // Definition of front cutout
  G4double           front_cutout_inner_rad;
  G4double           front_cutout_outer_rad;
  G4double           front_cutout_hx;
  
  // Define connector box cutout
  G4double            connector_thickness;
  G4double            connector_cutout_hx;
  G4double            connector_cutout_hy;
  G4double            connector_cutout_hz;
  
  // Define connector box
  G4double            connector_hx;
  G4double            connector_hy;
  G4double            connector_hz;
  
  G4Material*         fDefaultMaterial;
  G4double            fWorldSizeYZ;
  G4double            fWorldSizeX;

  G4Box*              fSolidWorld;
  G4LogicalVolume*    fLogicWorld;
  G4VPhysicalVolume*  fPhysiWorld;
  // Source Volume ---------------------------
  G4Material*         sourceContainerMaterial;
  G4double            fSourceVolume_x;
  G4double            fSourceVolume_y;
  G4double            fSourceVolume_z;
  G4double            fSourceVolumeDiameter;
  G4double            fSource_inner_rad;
  G4double            fSourceVolumeThickness;
  G4double            fSourceWindowThickness;
  G4Tubs*             fSourceVolumeTubs;
  G4LogicalVolume*    fLogicSource;
  G4VPhysicalVolume*  fPhysiSource;
  
  // New Cryo -----------------------
  G4Material*           cryoMaterial;
  G4Box*                cryoBox;
  G4LogicalVolume*      cryoLog;
  G4VPhysicalVolume*    cryoPhys;
  
  G4Material*           cryoCutoutMaterial;
  G4Box*                cryoCutoutBox;
  G4LogicalVolume*      cryoCutoutLog;
  G4VPhysicalVolume*    cryoCutoutPhys;
  
  G4Tubs*               frontCutoutTubs;
  G4SubtractionSolid*   cryoBoxSub_1;
  G4SubtractionSolid*   cryoBoxSub;
    
  
  G4Box*                connectorBox;
  G4LogicalVolume*      connectorLog;
  G4VPhysicalVolume*    connectorBoxPhys;
    
  G4Box*                connectorCutoutBox;
  G4LogicalVolume*      connectorCutoutLog;
  G4VPhysicalVolume*    connectorCutoutPhys;
    
  G4Box*                ifBox;
  G4LogicalVolume*      ifLog;
  G4VPhysicalVolume*    ifBoxPhys;
    
  G4Box*                ifCutoutBox;
  G4LogicalVolume*      ifCutoutLog;
  G4VPhysicalVolume*    ifCutoutPhys;
  //--------------------------------------
  G4Box*                fSolidCalor;
  G4LogicalVolume*      fLogicCalor;
  G4VPhysicalVolume*    fPhysiCalor;

  G4Box*                fSolidLayer;
  G4LogicalVolume*      fLogicLayer;
  G4VPhysicalVolume*    fPhysiLayer;

  G4Box*                fSolidAbsor[MaxAbsor];
  G4LogicalVolume*      fLogicAbsor[MaxAbsor];
  G4VPhysicalVolume*    fPhysiAbsor[MaxAbsor];
    
  G4double              CryoWallAbsorberStandOff;
  G4double              CalorOffset;
    
    //G4double		GeSizeYZ;
    //G4double		GeThickness;
    //G4Box*             solidGe;      //pointer to the solid Ge
    //G4LogicalVolume*   logicGe;      //pointer to the logical Ge
    //G4VPhysicalVolume* physiGe;      //pointer to the physical Ge

  //G4UniformMagField* fMagField;

  DetectorMessenger* fDetectorMessenger;

private:

  void DefineMaterials();
  void ComputeCalorParameters();
  G4VPhysicalVolume* ConstructCalorimeter();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

