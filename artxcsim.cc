#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "ActionInitialization.hh"
#include "G4UImanager.hh"

#include "DetectorConstruction.hh"
#include "XrayFluoPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"


#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"

#include "G4EmCalculator.hh"
#include "G4Gamma.hh"
#include "G4KleinNishinaCompton.hh"
#include "G4BetheHeitlerModel.hh"
#include "G4UnitsTable.hh"
#include <fstream>
#include <iostream>

int main(int argc,char** argv)
{
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }
  // Setting random	
  CLHEP::HepRandom::setTheSeed(time(0));//+getpid());
  // Construct the default run manager//
  // Try to use multitreaded regime with Nthreads  
    const G4int Nthreads = 1;
  #ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(Nthreads);
    printf("Multithreading");
  #else
    G4RunManager* runManager = new G4RunManager;
    printf("Single thread");
  #endif

  // Default detector geometry to be shared with other 
  // modules
  G4int sideN = 2; // Number of subpixels along detector axis
  G4double full_side = 28.56*mm; 
  G4double pix_side = full_side / sideN;
  G4double pix_depth = 1.*mm;
  // init detector
  G4VUserDetectorConstruction* detector = new DetectorConstruction(pix_side, pix_depth, sideN);
  runManager->SetUserInitialization(detector);
  // init physics
  G4VModularPhysicsList* physics = new PhysicsList;
  runManager->SetUserInitialization(physics);
  // set mandatory user action class
  //
  G4cout << "--->" <<G4endl;
  runManager->SetUserInitialization(new ActionInitialization);  

  // Initialize G4 kernel
  //
  runManager->Initialize();
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4cout << "--->" <<G4endl;
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  #ifdef G4VIS_USE
    delete visManager;
  #endif

  delete runManager;
}


