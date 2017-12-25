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
//
// $Id: SteppingAction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#ifndef G400Event_h
#define G400Event_h 1
#include "G400Event.hh"
#endif

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"


#include "G4VProcess.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(), eventaction(eventAction)					 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4Track* track         = aStep->GetTrack();
  G4StepPoint* prepoint = aStep->GetPreStepPoint();  
  G4StepPoint* postpoint = aStep->GetPostStepPoint();  
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4String vol_name = volume->GetName();
  if (vol_name.contains("CdTe") and edep > 0.0) 
    { 
      eventaction->AddStripe(vol_name, 0, edep);
      eventaction->AddTotEdep(edep);
    }  

  G4double ekin = track->GetKineticEnergy();	
  const G4DynamicParticle* step_particle = track->GetDynamicParticle();
  const G4ParticleDefinition* defin = step_particle->GetParticleDefinition();
  const G4String part_name = defin->GetParticleName();
  G4String ProcessName = "zzz";	
  if (track->GetCreatorProcess()==0 and ekin>1.*keV)
      {
        ProcessName = postpoint->GetProcessDefinedStep()->GetProcessName();
        if (ProcessName!="Transportation")
         {
           G4cout << part_name <<" " << ProcessName <<" " << ekin/keV << " " << G4endl;
         }
      }

    if( track->GetNextVolume() == 0 ) 
    {
      G4cout << "!leaving:"<< part_name <<" "<< ekin/keV << G4endl;  
      eventaction->AddLeavingParticle(part_name, ekin/keV);
    } 


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
