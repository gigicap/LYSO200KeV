
#include "G4ELINP_DetectorConstruction_NRSS.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4NistMessenger.hh"
#include "G4NistElementBuilder.hh"
#include "G4NistMaterialBuilder.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"

//constructor
G4ELINP_DetectorConstruction_NRSS::G4ELINP_DetectorConstruction_NRSS()
    : G4VUserDetectorConstruction()
{
}

G4ELINP_DetectorConstruction_NRSS::~G4ELINP_DetectorConstruction_NRSS()
{ }

G4VPhysicalVolume* G4ELINP_DetectorConstruction_NRSS::Construct()
{
    	
 	G4NistManager* nist = G4NistManager::Instance();
    	//materials    	
    	//Vacuum
    	G4double z = 7.;
    	G4double a = 14.007 * CLHEP::g/CLHEP::mole;
    	G4double density = CLHEP::universe_mean_density;
    	G4double pressure = 1.E-6 * 1.E-3 * CLHEP::bar;	//10-6 mbar
    	G4double temperature = 300. * CLHEP::kelvin;
    	G4Material *Vacuum = new G4Material("Vacuum", z, a, density, kStateGas, temperature, pressure);


//mother volume 
    G4Box* fNRSSEnvelopeSolid = new G4Box("NRSSEnvelopeSolid", 100*cm, 100*cm, 100*cm); 
    G4LogicalVolume *fNRSSEnvelopeLogic = new G4LogicalVolume(fNRSSEnvelopeSolid, Vacuum, "NRSSEnvelopeLogic");
    fNRSSEnvelopeLogic->SetVisAttributes(G4Colour(G4Colour(1.0,0.0,1.0)));  
    G4ThreeVector fNRSSEnvelopePositionVector = G4ThreeVector(0., 0., 0.);  
    G4VPhysicalVolume* fNRSSEnvelopePhysical = new G4PVPlacement(0, fNRSSEnvelopePositionVector, fNRSSEnvelopeLogic, "NRSSEnvelopePhysical", 0, false, 0); 
    
        


    	//assign matrials
       //G4double AC13 = 13.*g/mole;
       //G4double ZC13 = 6;
       //G4Element* C13  = new G4Element("Carbon13", "C", ZC13, AC13);
       //G4double densityC13 = 2.03*g/cm3;                                 
       //G4Material* target_mat = new G4Material("Target", densityC13, 1);
       G4Material* target_mat = nist->FindOrBuildMaterial("G4_Al");
       //target_mat->AddElement(C13, 1.);
    	
    
    	G4double fNRSSfrontScreenWidth = 5. * CLHEP::cm;
    	G4double fNRSSfrontScreenHeight = 20. * CLHEP::mm;
    	G4double fNRSSfrontScreenLength = 20. * CLHEP::mm;
  
  //first box 
      G4Box* fBoxSolid = new G4Box("BoxSolid", fNRSSfrontScreenLength*0.5, fNRSSfrontScreenHeight*0.5, fNRSSfrontScreenWidth*0.5);
      G4LogicalVolume *fBoxLogic = new G4LogicalVolume(fBoxSolid, target_mat, "BoxLogic");
      G4ThreeVector fBoxPositionVector = G4ThreeVector(0., 0., 0.);
      G4VPhysicalVolume* fBoxPhysical = new G4PVPlacement(0, fBoxPositionVector, fBoxLogic, "BoxPhysical", 0, false, 0); 
      new G4PVPlacement(0,         
                  fBoxPositionVector,    
                  fBoxLogic,             
                    "Box",                
                    fNRSSEnvelopeLogic,            
                    false,                  
                    0,                     
                    false);
 
  //second box (per vedere numeri ed energie in uscita dalla targhetta)
      G4double fNRSSfrontScreenWidth1 = 1. * CLHEP::mm;
      G4double fNRSSfrontScreenHeight1 = 20. * CLHEP::mm;
      G4double fNRSSfrontScreenLength1 = 20. * CLHEP::mm;


      G4Box* fBoxSolid1 = new G4Box("BoxSolid1", fNRSSfrontScreenLength1*0.5, fNRSSfrontScreenHeight1*0.5, fNRSSfrontScreenWidth1*0.5);
      G4LogicalVolume *fBoxLogic1 = new G4LogicalVolume(fBoxSolid1, Vacuum, "BoxLogic1");
      G4ThreeVector fBoxPositionVector1 = G4ThreeVector(0., 0., fNRSSfrontScreenWidth);
      G4VPhysicalVolume* fBoxPhysical1 = new G4PVPlacement(0, fBoxPositionVector1, fBoxLogic1, "BoxPhysical1", 0, false, 0); 
      new G4PVPlacement(0,         
                  fBoxPositionVector1,    
                  fBoxLogic1,             
                    "Box1",                
                    fNRSSEnvelopeLogic,            
                    false,                  
                    0,                     
                    false); 

   //LYSO 
  //central LYSO 
//LYSO
        G4Material* LYSO = new G4Material("LYSO", 7.1*g/cm3, 4);
        LYSO->AddElement(nist->FindOrBuildElement("Lu"), 71.45*perCent);
        LYSO->AddElement(nist->FindOrBuildElement("Si"), 6.37*perCent);
        LYSO->AddElement(nist->FindOrBuildElement("O"), 18.15*perCent);
        LYSO->AddElement(nist->FindOrBuildElement("Y"), 4.03*perCent);

        G4Material* LYSO_Ce = new G4Material("LYSO_Ce", 7.1*g/cm3, 2);
        LYSO_Ce->AddMaterial(LYSO, 99.81*perCent);
        LYSO_Ce->AddElement(nist->FindOrBuildElement("Ce"), 0.19*perCent);

      
        G4double centralSizex = 3 * cm;
        G4double centralSizey = 3 * cm;
        G4double centralDepth = 6 * cm;

        G4Box * centralSolid = new G4Box("centralSolid",centralSizex/2,centralSizey/2,centralDepth/2);

        G4LogicalVolume *centralLogic = new G4LogicalVolume(centralSolid, LYSO_Ce, "centralDetector");   
        
        G4VisAttributes* centralVisAttribute = new G4VisAttributes(G4Colour(0.,0.,1.));
        centralVisAttribute->SetForceSolid(true);
        centralLogic->SetVisAttributes(centralVisAttribute); 

        new G4PVPlacement(0,          
                    G4ThreeVector(0,0,50*cm +centralDepth/2),  
                    centralLogic,       
                    "CentralDetector",    
                    fNRSSEnvelopeLogic,            
                    false,                  
                    0,                     
                    false); 
   
    return fNRSSEnvelopePhysical;
}
