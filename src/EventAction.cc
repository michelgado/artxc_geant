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
// $Id: EventAction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "EventAction.hh"
#include "RunAction.hh"
#ifndef G400Event_h
#define G400Event_h 1
#include "G400Event.hh"
#endif

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction): G4UserEventAction(), runAct(runAction)
{
  printModulo = 1; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) 
    { 
      G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    }
  currentEvent.primary_name = "";
  currentEvent.primary_energy = 0.0;
  currentEvent.trigger = false;
  currentEvent.primx = -999;
  currentEvent.primy = -999;
  currentEvent.primz = -999;
  currentEvent.momx = -999;
  currentEvent.momy = -999;
  currentEvent.momz = -999;
  currentEvent.totalEdep = 0.0;
  currentEvent.stripes.clear();
  currentEvent.pads.clear();
  currentEvent.interactions.clear();
  currentEvent.leavingparts.clear();


  G4double primary_kin = evt->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();
  G4String primary_name = evt->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetParticleName();
  currentEvent.primary_name = primary_name;
  currentEvent.primary_energy = primary_kin/MeV;
  currentEvent.primx = evt->GetPrimaryVertex()->GetX0();
  currentEvent.primy = evt->GetPrimaryVertex()->GetY0();
  currentEvent.primz = evt->GetPrimaryVertex()->GetZ0();
  currentEvent.momx = evt->GetPrimaryVertex()->GetPrimary()->GetPx();
  currentEvent.momy = evt->GetPrimaryVertex()->GetPrimary()->GetPy();
  currentEvent.momz = evt->GetPrimaryVertex()->GetPrimary()->GetPz();
  G4cout << "---> Primary particle "<< primary_name <<" with energy, keV: " << primary_kin/keV << G4endl;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();
  G4cout << "--->Total energy deposit:" << currentEvent.totalEdep/keV << G4endl;	
  G4cout << "---> End of event: " << evtNb << G4endl;	
  while ((int)currentEvent.leavingparts.size() > 0)
          {
            G400phEvent temp = currentEvent.leavingparts.back();
            G4cout <<" lp_"<<temp.procname << " " << temp.ekin << std::endl;
            currentEvent.leavingparts.pop_back();
          }
  runAct->fillPerEvent(currentEvent);
}  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::AddStripe(G4String smothername , G4int snumber, G4double sedep)
{ 
  G400workedStripe inserting, temp;
  inserting.motherName = smothername;
  inserting.number = snumber;  
  inserting.edep = sedep/MeV;
  if ((int)currentEvent.stripes.size() == 0) 
      currentEvent.stripes.push_back(inserting);
  else
    {
      temp  = currentEvent.stripes.back();
      currentEvent.stripes.pop_back();
      if (temp.motherName == inserting.motherName and temp.number == inserting.number)
        {  
          sedep = (sedep/MeV + temp.edep);
          inserting.edep = sedep;
          currentEvent.stripes.push_back(inserting);          
        }   
      else 
        {
          currentEvent.stripes.push_back(temp);            
          G4bool found = false;  
          for(int y=0; y<currentEvent.stripes.size(); y++) 
            {
              temp  = currentEvent.stripes[y];
              if (temp.motherName == inserting.motherName and temp.number == inserting.number)
                { 
                  currentEvent.stripes.erase(currentEvent.stripes.begin()+y);
                  sedep = (sedep/MeV + temp.edep);
                  inserting.edep = sedep;
                  currentEvent.stripes.push_back(inserting);          
                  found = true;
                  break;
                }   
            }  
          if (!found)
            {  
              currentEvent.stripes.push_back(inserting);
            }  
        }  
    }    
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::AddTotEdep(G4double sedep)
{
  currentEvent.totalEdep+=sedep;
}

void  EventAction::AddLeavingParticle(G4String pname, G4double partekin)
{
  G400phEvent inserting;
  inserting.procname = pname;
  inserting.ekin = partekin;
  currentEvent.leavingparts.push_back(inserting);
}


