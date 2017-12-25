//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "globals.hh"
#include <vector>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct G400particle
    {
    G4String name;
    G4double energy;
    G4double time;
    };
struct G400Edep
    {
    G4double energy;
    G4double time;
    G4String pname;
    };
    
struct G400workedStripe
    {
    G4String motherName;
    G4int number;
    G4double edep;
    };

struct G400workedPad
    {
    G4String name;
    G4int numberY;
    G4int numberZ;
    G4double edep;
    };

struct G400phEvent
    {
    G4String procname;
    G4double ekin;
    };

    
struct G400Event{
  G4String   primary_name; //Name of a primary particle
  G4double primary_energy; //Energy of a primary particle
  std::vector<G400workedStripe> stripes;   
  std::vector<G400workedPad> pads;   
  std::vector<G400phEvent> interactions;   
  std::vector<G400phEvent> leavingparts;  
  G4double S1time, S2time;
  G4double primx, primy, primz, momx, momy, momz;
  G4bool trigger;
  G4double totalEdep;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

