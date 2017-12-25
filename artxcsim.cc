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
  #ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(20);
    printf("Parallel!");
  #else
    G4RunManager* runManager = new G4RunManager;
    printf("Single-treaded!");
  #endif

  // set mandatory initialization classes//
  G4double full_side = 28.56*mm;
  G4double pix_side = full_side / 96.;
  G4double pix_depth = 1.*mm;
  G4double electrode_layer_depth = 0.05*mm/2; //50 micrimeters
  G4double stripe_width = 0.52*mm;
  G4double stripe_spacing = 0.075*mm; 
  G4int sideN = 96;
  G4cout << "--->" <<G4endl;
  G4VUserDetectorConstruction* detector = new DetectorConstruction(pix_side, pix_depth, sideN);
  runManager->SetUserInitialization(detector);
  //
  G4cout << "--->" <<G4endl;
  G4VModularPhysicsList* physics = new PhysicsList;
  runManager->SetUserInitialization(physics);

  // set mandatory user action class
  //
  G4cout << "--->" <<G4endl;
  runManager->SetUserInitialization(new ActionInitialization);  
  //G4VUserPrimaryGeneratorAction* gen_action = new PrimaryGeneratorAction(pix_side);
  //runManager->SetUserAction(gen_action);
  //runManager->SetUserAction(new RunAction);
  //runManager->SetUserAction(new EventAction);
  //runManager->SetUserAction(new SteppingAction);

  // Initialize G4 kernel
  //
  runManager->Initialize();
  G4cout << "rm--->" <<G4endl;
// Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  // Get the pointer to the UI manager and set verbosities
  //
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
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  }
  #ifdef G4VIS_USE
    delete visManager;
  #endif
  delete runManager;
  return 0;
}


