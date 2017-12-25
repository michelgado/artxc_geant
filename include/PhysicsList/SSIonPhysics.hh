//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Title:  Ion Physics for Space Sciences Physics List                    //
//  Date:   22 March 2005                                                  //
//  Author: D.H. Wright (SLAC)                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

#ifndef SSIonPhysics_h
#define SSIonPhysics_h 1

#include "G4VPhysicsConstructor.hh"


class SSIonPhysics : public G4VPhysicsConstructor
{
  public: 
    SSIonPhysics(const G4String& name ="ion");
    virtual ~SSIonPhysics();

    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

};

#endif





