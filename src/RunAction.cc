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
// $Id: RunAction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "RunAction.hh"
#ifndef G400Event_h
#define G400Event_h 1
#include "G400Event.hh"
#endif


#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <iostream>
#include <algorithm>    // std::reverse
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
  {
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent(G400Event current)
{
  event_store.push_back(current);
  if ( (int)event_store.size() > 50)
  {
    std::ofstream myfile;
    myfile.open ("events.eventdata",  std::ios::out | std::ios::app);
    while ((int)event_store.size() > 0)
      { 
        G400Event to_wr = event_store.back();
        event_store.pop_back();
        if ((int)to_wr.stripes.size() == 0)
			{continue;}
        myfile << to_wr.primary_name << " " << to_wr.primary_energy << " " << to_wr.primy <<" "<<to_wr.primz;
	while ((int)to_wr.stripes.size() > 0)
          {
            G400workedStripe temp = to_wr.stripes.back();
            myfile << " " << temp.motherName << " "<<temp.edep;
            to_wr.stripes.pop_back();
          }

        while ((int)to_wr.leavingparts.size() > 0)
          {
            G400phEvent temp = to_wr.leavingparts.back();
            myfile <<" lp_"<<temp.procname << " " << temp.ekin << std::endl;
            to_wr.leavingparts.pop_back();
          }

	myfile << std::endl;

      }
    myfile.close();
  }      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) 
      return;
  if ( (int)event_store.size() > 0)
    {
      std::ofstream myfile;
      myfile.open ("events.eventdata",  std::ios::out | std::ios::app);
      while ((int)event_store.size() > 0)
        {
          G400Event to_wr = event_store.back();
          event_store.pop_back();
		  if ((int)to_wr.stripes.size() == 0)
			{continue;}
        myfile << to_wr.primary_name << " " << to_wr.primary_energy << " " << to_wr.primy <<" "<<to_wr.primz;
	while ((int)to_wr.stripes.size() > 0)
          {
            G400workedStripe temp = to_wr.stripes.back();
            myfile << " " << temp.motherName << " "<<temp.edep;
            to_wr.stripes.pop_back();
          }
        while ((int)to_wr.leavingparts.size() > 0)
          {
            G400phEvent temp = to_wr.leavingparts.back();
            myfile <<" lp_"<<temp.procname << " " << temp.ekin << std::endl;
            to_wr.leavingparts.pop_back();
          }

	myfile << std::endl;
        }
      myfile.close();
   }      
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
