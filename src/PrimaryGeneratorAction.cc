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
// $Id: ExN01PrimaryGeneratorAction.cc,v 1.6 2006-06-29 17:47:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "PrimaryGeneratorAction.hh"

#include "Randomize.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
  particleGun->SetParticlePosition(G4ThreeVector(-1.0*mm, 0.0, 0.0));
  particleGun->SetParticleEnergy(60.*keV);  
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double start_x = -45.25*mm;
  G4double start_y = ((G4UniformRand()*2-1)*0.5*mm);
  G4double start_z = ((G4UniformRand()*2-1)*0.5*mm);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
  G4double cosTheta = -1.0+2.0*G4UniformRand();
  G4double phi = twopi * G4UniformRand();
  G4double sinTheta = sqrt(1-cosTheta*cosTheta);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(sinTheta * cos(phi) , sinTheta * sin(phi), cosTheta));

  particleGun->SetParticlePosition(G4ThreeVector(start_x,start_y,start_z));
  particleGun->GeneratePrimaryVertex(anEvent);
}


