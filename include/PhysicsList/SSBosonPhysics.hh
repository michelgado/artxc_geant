//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Title:  Boson Physics for Space Sciences Phyiscs List                  //
//  Date:   22 March 2005                                                  //
//  Author: D.H. Wright (SLAC)                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

#ifndef SSBosonPhysics_h
#define SSBosonPhysics_h 1

#include "G4VPhysicsConstructor.hh"

class SSBosonPhysics : public G4VPhysicsConstructor
{
  public: 
    SSBosonPhysics(const G4String& name = "boson");
    virtual ~SSBosonPhysics();

    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

};

#endif








