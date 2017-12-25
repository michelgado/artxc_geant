//
////////////////////////////////////////////////////////////////
//                                                            //
//  Title:  Decay Physics for Space Sciences Physics List     //
//  Date:   20 March 2005                                     //
//  Author: D.H. Wright (SLAC)                                //
//                                                            //
////////////////////////////////////////////////////////////////
//

#ifndef SSDecayPhysics_h
#define SSDecayPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "G4Decay.hh"


class SSDecayPhysics : public G4VPhysicsConstructor
{
  public: 
    SSDecayPhysics(const G4String& name = "decay");
    virtual ~SSDecayPhysics();

    // This method will be invoked in the Construct() method.
    // each particle type will be instantiated
    virtual void ConstructParticle();

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

  protected:
    G4Decay fDecayProcess;
};


#endif








