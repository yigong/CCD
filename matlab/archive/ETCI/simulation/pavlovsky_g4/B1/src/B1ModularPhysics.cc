#include "B1ModularPhysics.hh"

#include "G4RegionStore.hh"

#include "G4ProcessManager.hh"
#include "G4UserLimits.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4PhysicsListHelper.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"
#include "G4EmProcessOptions.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4DecayPhysics.hh"

#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4RadioactiveDecayPhysics.hh"

B1ModularPhysics::B1ModularPhysics() : G4VModularPhysicsList()
{
//    G4LossTableManager::Instance();
    defaultCutValue = 10.0*nanometer;
    SetVerboseLevel(1);

    emPhysicsList = new G4EmLivermorePhysics();

    decPhysicsList = new G4RadioactiveDecayPhysics();

    particleList = new G4DecayPhysics();
}

B1ModularPhysics::~B1ModularPhysics()
{
}

void B1ModularPhysics::ConstructParticle()
{
    particleList->ConstructParticle();
}

void B1ModularPhysics::ConstructProcess()
{

    // transportation
    //
    AddTransportation();

    // electromagnetic physics list
    //
    emPhysicsList->ConstructProcess();

    // particle physics list
    //
    particleList->ConstructProcess();

    // decay physics list
    //
    decPhysicsList->ConstructProcess();


}

void B1ModularPhysics::SetCuts()
{

    if (verboseLevel >0) {
      G4cout << "PhysicsList::SetCuts:";
      G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
    }

    // set cut values for gamma at first and for e- second and next for e+,
    // because some processes for e+/e- need cut values for gamma
    SetCutValue(defaultCutValue, "gamma");
    SetCutValue(defaultCutValue, "e-");
    SetCutValue(defaultCutValue, "e+");

    G4ProductionCuts *opt_cut = new G4ProductionCuts();
    opt_cut->SetProductionCut(defaultCutValue/100.,"gamma");
    opt_cut->SetProductionCut(defaultCutValue/100.,"e-");
    opt_cut->SetProductionCut(defaultCutValue/100.,"e+");
    opt_cut->SetProductionCut(defaultCutValue/100.,"proton");


    G4Region* ccdLogReg = (G4RegionStore::GetInstance())->GetRegion("Ccd");
    ccdLogReg->SetProductionCuts(opt_cut);

    G4Region *ccdMountBoardLogReg = (G4RegionStore::GetInstance())->GetRegion("ccdMountBoard");
    ccdMountBoardLogReg->SetProductionCuts(opt_cut);

    G4Region *alNitrideBackLogReg = (G4RegionStore::GetInstance())->GetRegion("alNitrideBack");
    alNitrideBackLogReg->SetProductionCuts(opt_cut);

    G4Region *pfbLogReg = (G4RegionStore::GetInstance())->GetRegion("pfb");
    pfbLogReg->SetProductionCuts(opt_cut);

    G4Region *coldPlateLogReg = (G4RegionStore::GetInstance())->GetRegion("coldPlate");
    coldPlateLogReg->SetProductionCuts(opt_cut);

    G4Region *copperBarLogReg = (G4RegionStore::GetInstance())->GetRegion("copperBar");
    copperBarLogReg->SetProductionCuts(opt_cut);

    G4Region *cryoLogReg = (G4RegionStore::GetInstance())->GetRegion("cryoLog");
    cryoLogReg->SetProductionCuts(opt_cut);

    G4EmProcessOptions emOptions;

    emOptions.SetFluo(true); // To activate deexcitation processes and fluorescence
    emOptions.SetAuger(true); // To activate Auger effect if deexcitation is activated
    emOptions.SetPIXE(true); // To activate Particle Induced X-Ray Emission (PIXE)

    emOptions.SetDeexcitationActiveRegion("Ccd");

    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250.*eV, 50.*MeV);

    if (verboseLevel>0) DumpCutValuesTable();

}

//void B1ModularPhysics::ConstructDecay()
//{

//}

//void B1ModularPhysics::ConstructEM()
//{

//}
