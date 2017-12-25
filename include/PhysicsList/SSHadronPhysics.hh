//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Title:  Hadron Physics for Space Sciences Physics List                 //
//  Date:   20 March 2005                                                  //
//  Author: D.H. Wright (SLAC)                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

#ifndef SSHadronPhysics_h
#define SSHadronPhysics_h 1

#include "G4VPhysicsConstructor.hh"

class SSHadronPhysics : public G4VPhysicsConstructor
{
  public: 
    SSHadronPhysics(const G4String& name ="hadron");
    virtual ~SSHadronPhysics();

    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

};

#endif





