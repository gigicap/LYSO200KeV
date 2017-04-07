#ifndef G4ELINP_DetectorConstruction_NRSS_h
#define G4ELINP_DetectorConstruction_NRSS_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"



#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4NistMessenger.hh"
#include "G4NistElementBuilder.hh"
#include "G4NistMaterialBuilder.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"



class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class G4ELINP_DetectorConstruction_NRSS : public G4VUserDetectorConstruction
{
  public:    
    G4ELINP_DetectorConstruction_NRSS();
    virtual ~G4ELINP_DetectorConstruction_NRSS();

    virtual G4VPhysicalVolume* Construct();  
    

};

#endif
   
