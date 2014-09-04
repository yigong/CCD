#include "B1PhysicsList.hh"

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

#include "G4Decay.hh"

B1PhysicsList::B1PhysicsList()
{
    defaultCutValue = 100.0*nanometer;
    SetVerboseLevel(1);
}

B1PhysicsList::~B1PhysicsList()
{
}

void B1PhysicsList::ConstructParticle()
{
    // In this method, static member functions should be called
    // for all particles which you want to use.
    // This ensures that objects of these particle types will be
    // created in the program.

    G4BosonConstructor  pBosonConstructor;
    pBosonConstructor.ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();

    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();
}

void B1PhysicsList::ConstructProcess()
{
    AddTransportation();
    ConstructEM();
    ConstructDecay();
}

void B1PhysicsList::SetCuts()
{
    if (verboseLevel >0){
      G4cout << "PhysicsList::SetCuts:";
      G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
    }

    // set cut values for gamma at first and for e- second and next for e+,
    // because some processes for e+/e- need cut values for gamma
    //
    SetCutValue(defaultCutValue, "gamma");
    SetCutValue(defaultCutValue, "e-");
    SetCutValue(defaultCutValue, "e+");
    SetCutValue(defaultCutValue, "proton");

    G4ProductionCuts *opt_cut = new G4ProductionCuts();
//    fRadiatorCuts->SetProductionCut(fGammaCut, idxG4GammaCut);
//    fRadiatorCuts->SetProductionCut(fElectronCut, idxG4ElectronCut);
//    fRadiatorCuts->SetProductionCut(fPositronCut, idxG4PositronCut);

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

    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100.*eV, 1.*GeV);

    G4cout << "CUTS SET";

    if (verboseLevel>0) DumpCutValuesTable();
}

void B1PhysicsList::ConstructDecay()
{
//    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

//    // Add Decay Process
//    G4Decay* theDecayProcess = new G4Decay();
//    theParticleIterator->reset();
//    while( (*theParticleIterator)() ){
//      G4ParticleDefinition* particle = theParticleIterator->value();
//      if (theDecayProcess->IsApplicable(*particle)) {
//        ph->RegisterProcess(theDecayProcess, particle);
//      }
//    }
    G4Decay* theDecayProcess = new G4Decay();
    theParticleIterator->reset();

    while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if (theDecayProcess->IsApplicable(*particle))
      {
        pmanager ->AddProcess(theDecayProcess);

        // set ordering for PostStepDoIt and AtRestDoIt

        pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
        pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
      }
    }
}

#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4UniversalFluctuation.hh"

void B1PhysicsList::ConstructEM()
{

    // Add standard EM Processes

    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      //Applicability range for Livermore models
      //for higher energies, the Standard models are used
      G4double highEnergyLimit = 1*GeV;

      if (particleName == "gamma") {
        // gamma

        G4PhotoElectricEffect* phot = new G4PhotoElectricEffect();
        G4LivermorePhotoElectricModel*
        photModel = new G4LivermorePhotoElectricModel();
        photModel->SetHighEnergyLimit(highEnergyLimit);
        phot->AddEmModel(0, photModel);
        pmanager->AddDiscreteProcess(phot);

        G4ComptonScattering* compt = new G4ComptonScattering();
        G4LivermoreComptonModel*
        comptModel = new G4LivermoreComptonModel();
        comptModel->SetHighEnergyLimit(highEnergyLimit);
        compt->AddEmModel(0, comptModel);
        pmanager->AddDiscreteProcess(compt);

        G4GammaConversion* conv = new G4GammaConversion();
        G4LivermoreGammaConversionModel*
        convModel = new G4LivermoreGammaConversionModel();
        convModel->SetHighEnergyLimit(highEnergyLimit);
        conv->AddEmModel(0, convModel);
        pmanager->AddDiscreteProcess(conv);

        G4RayleighScattering* rayl = new G4RayleighScattering();
        G4LivermoreRayleighModel*
        raylModel = new G4LivermoreRayleighModel();
        raylModel->SetHighEnergyLimit(highEnergyLimit);
        rayl->AddEmModel(0, raylModel);
        pmanager->AddDiscreteProcess(rayl);

      } else if (particleName == "e-") {
        //electron

        G4eIonisation* eIoni = new G4eIonisation();
        G4LivermoreIonisationModel*
        eIoniModel = new G4LivermoreIonisationModel();
        eIoniModel->SetHighEnergyLimit(highEnergyLimit);
        eIoni->AddEmModel(0, eIoniModel, new G4UniversalFluctuation() );
        pmanager->AddProcess(eIoni,                   -1,-1, 1);

        G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
        G4LivermoreBremsstrahlungModel*
        eBremModel = new G4LivermoreBremsstrahlungModel();
        eBremModel->SetHighEnergyLimit(highEnergyLimit);
        eBrem->AddEmModel(0, eBremModel);
        pmanager->AddProcess(eBrem,                   -1,-1, 2);

      } else if (particleName == "e+") {
        //positron
        pmanager->AddProcess(new G4eIonisation,       -1,-1, 1);
        pmanager->AddProcess(new G4eBremsstrahlung,   -1,-1, 2);
        pmanager->AddProcess(new G4eplusAnnihilation,  0,-1, 3);

      } else if( particleName == "mu+" ||
                 particleName == "mu-"    ) {
        //muon
        pmanager->AddProcess(new G4MuIonisation,      -1,-1, 1);
        pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1, 2);
        pmanager->AddProcess(new G4MuPairProduction,  -1,-1, 3);

      } else if( particleName == "alpha" || particleName == "GenericIon" ) {
        pmanager->AddProcess(new G4ionIonisation,     -1,-1, 1);

      } else if ((!particle->IsShortLived()) &&
             (particle->GetPDGCharge() != 0.0) &&
             (particle->GetParticleName() != "chargedgeantino")) {
        //all others charged particles except geantino
        pmanager->AddProcess(new G4hIonisation,       -1,-1, 1);
      }
    }
//    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

//    theParticleIterator->reset();
//    while( (*theParticleIterator)() ){
//      G4ParticleDefinition* particle = theParticleIterator->value();
//      G4String particleName = particle->GetParticleName();

//      if (particleName == "gamma") {
//        // gamma
//        ph->RegisterProcess(new G4PhotoElectricEffect, particle);
//        ph->RegisterProcess(new G4ComptonScattering,   particle);
//        ph->RegisterProcess(new G4GammaConversion,     particle);

//      } else if (particleName == "e-") {
//        //electron
//        ph->RegisterProcess(new G4eMultipleScattering, particle);
//        ph->RegisterProcess(new G4eIonisation,         particle);
//        ph->RegisterProcess(new G4eBremsstrahlung,     particle);

//      } else if (particleName == "e+") {
//        //positron
//        ph->RegisterProcess(new G4eMultipleScattering, particle);
//        ph->RegisterProcess(new G4eIonisation,         particle);
//        ph->RegisterProcess(new G4eBremsstrahlung,     particle);
//        ph->RegisterProcess(new G4eplusAnnihilation,   particle);

//      } else if( particleName == "mu+" ||
//                 particleName == "mu-"    ) {
//        //muon
//        ph->RegisterProcess(new G4MuMultipleScattering, particle);
//        ph->RegisterProcess(new G4MuIonisation,         particle);
//        ph->RegisterProcess(new G4MuBremsstrahlung,     particle);
//        ph->RegisterProcess(new G4MuPairProduction,     particle);

//      } else if( particleName == "proton" ||
//                 particleName == "pi-" ||
//                 particleName == "pi+"    ) {
//        //proton
//        ph->RegisterProcess(new G4hMultipleScattering, particle);
//        ph->RegisterProcess(new G4hIonisation,         particle);
//        ph->RegisterProcess(new G4hBremsstrahlung,     particle);
//        ph->RegisterProcess(new G4hPairProduction,     particle);

//      } else if( particleName == "alpha" ||
//             particleName == "He3" )     {
//        //alpha
//        ph->RegisterProcess(new G4hMultipleScattering, particle);
//        ph->RegisterProcess(new G4ionIonisation,       particle);

//      } else if( particleName == "GenericIon" ) {
//        //Ions
//        ph->RegisterProcess(new G4hMultipleScattering, particle);
//        ph->RegisterProcess(new G4ionIonisation,       particle);

//        } else if ((!particle->IsShortLived()) &&
//             (particle->GetPDGCharge() != 0.0) &&
//             (particle->GetParticleName() != "chargedgeantino")) {
//        //all others charged particles except geantino
//        ph->RegisterProcess(new G4hMultipleScattering, particle);
//        ph->RegisterProcess(new G4hIonisation,         particle);
//      }
//    }
    G4EmProcessOptions emOptions;

    emOptions.SetFluo(true); // To activate deexcitation processes and fluorescence
    emOptions.SetAuger(true); // To activate Auger effect if deexcitation is activated
    emOptions.SetPIXE(true); // To activate Particle Induced X-Ray Emission (PIXE)

    emOptions.SetDeexcitationActiveRegion("Ccd");
}
