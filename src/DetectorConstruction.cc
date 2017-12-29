#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Polycone.hh"
#include "G4RotationMatrix.hh"

#include "G4PhysicalConstants.hh"
#include "G4NistManager.hh"

#include "globals.hh"
#include <sstream>

DetectorConstruction::DetectorConstruction(G4double pixel_side, G4double pixel_depth, G4int pixel_side_N)
 :  space_log(0), space_phys(0)
{
        pix_side=pixel_side;
        pix_depth=pixel_depth;
        pix_side_num=pixel_side_N;	
}
//...oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo...//  
DetectorConstruction::~DetectorConstruction()
{
}
//...oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo......oooOOO000OOOooo...//
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Init elements and materials for detector	
  G4Element* elemCd   = new G4Element( "Cadmium"  , "Cd"    , 48 , 112.414    * g/mole );
  G4Element* elemTe   = new G4Element( "Tellurium"  , "Te"    , 52 , 127.6    * g/mole );
  G4Element* elemNi   = new G4Element( "Nickel"  , "Ni"    , 28 , 58.6934   * g/mole );
  G4Element* elemCu   = new G4Element( "Copper"  , "Cu"    , 29 , 63.546   * g/mole );


  G4Material* materVacuum = new G4Material( "Vacuum" ,  1 , 1.01*g/mole, universe_mean_density , kStateGas , 2.73*kelvin , 3.e-18*pascal );
  G4Material* materCdTe   = new G4Material( "CdTe     "            , 5.85*g/cm3 , 2     );
  G4Material* materMonel   = new G4Material( "Monel"            , 8.8*g/cm3 , 2     );

  materCdTe   -> AddElement( elemCd  , 1);
  materCdTe   -> AddElement( elemTe  , 1);  

  materMonel   -> AddElement( elemNi  , 2);  
  materMonel   -> AddElement( elemCu  , 1);  


  G4NistManager* man = G4NistManager::Instance();
  G4Material* materAu  = man->FindOrBuildMaterial("G4_Au");
  G4Material* materPt  = man->FindOrBuildMaterial("G4_Pt");  
  G4Material* materTi  = man->FindOrBuildMaterial("G4_Ti");  
  G4Material* materAl  = man->FindOrBuildMaterial("G4_Al");  
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //Defining volume of world//
  G4double space_x = 1.0*m;
  G4double space_y = 1.0*m;
  G4double space_z = 1.0*m;
  G4Box* space_box = new G4Box("space_box", space_x,space_y,space_z);
  space_log = new G4LogicalVolume(space_box, materVacuum, "space_log",0,0,0);
  space_phys = new G4PVPlacement(0,G4ThreeVector(),space_log,"space",0,false,checkOverlaps);

  G4double zs[] = {0, 4.5*mm, 4.5*mm, 5.0*mm};
  G4double rins[] = {0.,0.,0.,1.*mm};
  G4double routs[] = {4.*mm,4.*mm,4.*mm,4.*mm};
  G4Polycone* radsource = new G4Polycone("Radsource", 0., 2*pi, 4, zs, rins,routs);
  G4LogicalVolume* radsource_log = new G4LogicalVolume(radsource, materMonel, "Radsource", 0, 0 ,0);
  G4RotationMatrix * rotY = new G4RotationMatrix;
  rotY->rotateY( 270.0 * degree );
  new G4PVPlacement(rotY, G4ThreeVector(-5.0*cm, 0., 0.), radsource_log, "Radsource", space_log, 0, 0);

  //General sizes
  G4double electrode_layer_depth = (500*1e-6)*mm/2; //500 nanometers
  G4int side_size = pix_side_num; //number of pixels in detector side
  G4int thNlayers = 1;
  G4double thCdTe_depth = (0.05*pix_depth/2.)/thNlayers;    
  G4double CdTe_side = pix_side/2.;
  G4double full_side = 28.56*mm;

  //All detector lie in positive x's
  G4double curr_x = 0.;
  
  //Au-Pt layers on top of detector
  G4Box* electrode_box = new G4Box("electrode_box", space_x,space_y,space_z);
  G4Box* electrode_layer_box = new G4Box("electrode_layer_box", electrode_layer_depth, full_side/2., full_side/2.);
  G4LogicalVolume* au_electrode_layer_log = new G4LogicalVolume(electrode_layer_box, materAu, "Au_layer", 0, 0 ,0);
  G4LogicalVolume* pt_electrode_layer_log = new G4LogicalVolume(electrode_layer_box, materPt, "Pt_layer", 0, 0 ,0);
  new G4PVPlacement(0, G4ThreeVector(curr_x+electrode_layer_depth, 0.,
                       0.), au_electrode_layer_log, "Au_anode", space_log, 0, 0);
  curr_x+=2.*electrode_layer_depth;
  new G4PVPlacement(0, G4ThreeVector(curr_x+electrode_layer_depth, 0.,
                       0.), pt_electrode_layer_log, "Pt_anode", space_log, 0, 0);
  curr_x+=2.*electrode_layer_depth;


  //ten layers of thin CdTe pixel
  G4Box* thCdTe_pixel_box = new G4Box("CdTe", thCdTe_depth, CdTe_side, CdTe_side);
  G4LogicalVolume* thCdTe_pixel_log = new G4LogicalVolume(thCdTe_pixel_box, materCdTe, "CdTe_pixel", 0, 0 ,0);

  //Create pixel array        
  G4String curr_name;
  for (int n=0; n<thNlayers; n++) //X-axis iterator
  {
  G4cout << "--->th" << curr_x << G4endl;
  for (int ii=0; ii<side_size; ii++) //Y-axis iterator
    {
      for (int kk=0; kk<side_size; kk++) //Z-axis iterator  
            {
              //Generate number of element  
              std::stringstream ss;
              ss<<(ii*side_size+kk)<<"_"<<n;
              ss>>curr_name;
              curr_name = G4String("thCdTe_")+curr_name;
              // Place pixel     
              new G4PVPlacement(0, G4ThreeVector(curr_x+thCdTe_depth, (2*CdTe_side*(ii - side_size*0.5 + 0.5)),
                               (2*CdTe_side*(kk - side_size*0.5 + 0.5))), thCdTe_pixel_log, curr_name, space_log, 0, 0);
            }
        }
  curr_x+=2.*thCdTe_depth;
  }

  //Now placing regular thick pixels
  G4int Nlayers = 5;
  G4double CdTe_depth = (0.95*pix_depth/2.)/Nlayers;    


  //CdTe pixel
  G4Box* CdTe_pixel_box = new G4Box("CdTe", CdTe_depth, CdTe_side, CdTe_side);
  G4LogicalVolume* CdTe_pixel_log = new G4LogicalVolume(CdTe_pixel_box, materCdTe, "CdTe_pixel", 0, 0 ,0);

  //Create pixel array        

  for (int n=0; n<Nlayers; n++) //X-axis iterator
  {
  G4cout << "--->tk" << curr_x << G4endl;
  for (int ii=0; ii<side_size; ii++) //Y-axis iterator
    {
      for (int kk=0; kk<side_size; kk++) //Z-axis iterator  
            {
              //Generate number of element  
              std::stringstream ss;
              ss<<(ii*side_size+kk)<<"_"<<n;
              ss>>curr_name;
              curr_name = G4String("CdTe_")+curr_name;
              // Place pixel     
              new G4PVPlacement(0, G4ThreeVector(curr_x+CdTe_depth, (2*CdTe_side*(ii - side_size*0.5 + 0.5)),
                               (2*CdTe_side*(kk - side_size*0.5 + 0.5))), CdTe_pixel_log, curr_name, space_log, 0, 0);
//	      G4cout << "--->" <<(2*CdTe_side*(ii - side_size*0.5 + 0.5))<< " " << (2*CdTe_side*(kk - side_size*0.5 + 0.5)) << G4endl;
            }
        }
  curr_x+=2.*CdTe_depth;
  }



  G4LogicalVolume* al_electrode_layer_log = new G4LogicalVolume(electrode_layer_box, materAl, "Al_layer", 0, 0 ,0);
  G4LogicalVolume* ti_electrode_layer_log = new G4LogicalVolume(electrode_layer_box, materTi, "Ti_layer", 0, 0 ,0);
  new G4PVPlacement(0, G4ThreeVector(curr_x+electrode_layer_depth, 0.,
                       0.), al_electrode_layer_log, "Al_cathode", space_log, 0, 0);
  curr_x+=2.*electrode_layer_depth;
  new G4PVPlacement(0, G4ThreeVector(curr_x+electrode_layer_depth, 0.,
                       0.), ti_electrode_layer_log, "Ti_cathode", space_log, 0, 0);
  curr_x+=2.*electrode_layer_depth;
  new G4PVPlacement(0, G4ThreeVector(curr_x+electrode_layer_depth, 0.,
                       0.), au_electrode_layer_log, "Au_cathode", space_log, 0, 0);
  curr_x+=2.*electrode_layer_depth;



  //Set up visualisation 
  space_log->SetVisAttributes (G4VisAttributes::Invisible);
  
  G4VisAttributes* CdTe_vis= new G4VisAttributes(G4Colour(1.0,0.5,0.0));
  G4VisAttributes* rs_vis= new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  CdTe_vis->SetVisibility(true);
  CdTe_pixel_log->SetVisAttributes(CdTe_vis);
  radsource_log->SetVisAttributes(rs_vis);
  //
  //always return the physical World
  //
  return space_phys;
}
