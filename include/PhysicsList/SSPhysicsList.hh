//
////////////////////////////////////////////////////////////////
//                                                            //
//  Title:  Space Sciences Physics List                       //
//  Date:   20 March 2005                                     //
//  Author: D.H. Wright (SLAC)                                //
//                                                            //
////////////////////////////////////////////////////////////////
//

#ifndef SSPhysicsList_h
#define SSPhysicsList_h 1

#include "G4VModularPhysicsList.hh"


class SSPhysicsList: public G4VModularPhysicsList
{
public:
  SSPhysicsList();
  virtual ~SSPhysicsList();
  
  virtual void SetCuts();
};

#endif



