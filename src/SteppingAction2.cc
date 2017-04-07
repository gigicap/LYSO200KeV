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
/// \file electromagnetic/TestEm4/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
// $Id: SteppingAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "G4ELINP_DetectorConstruction_NRSS.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4OpticalPhoton.hh"
#include "G4UImanager.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* EvAct, PrimaryGeneratorAction * primGen)
    :G4UserSteppingAction(),fEventAction(EvAct),fPrimaryGeneratorAction(primGen)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "Analysis.hh"
#include <cmath>
void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    using namespace CLHEP;
    using namespace std;
    G4String CPName;
 
   G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
 
    G4VPhysicalVolume * volume = aStep->GetPreStepPoint()->GetPhysicalVolume();

    //cout<<volume->GetName()<<endl;

    if(volume->GetName() == "Box1"){
        //cout<<"Got hit"<<endl;
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
			if(aStep->GetPreStepPoint()->GetStepStatus()==fGeomBoundary){		//se entra nel rivelatore
                //fEventAction->addDetectorNh(0,1);
				G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();					//energia cinetica incidente
                G4ThreeVector direction = aStep->GetPreStepPoint()->GetMomentumDirection();
                G4double time = aStep->GetPreStepPoint()->GetGlobalTime();                      //tempo 
                G4int PDGcode = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();          //PDG code
                //traccia dello step
                G4int CPcode = 0;   
                G4int Nstep = aStep->GetTrack()->GetCurrentStepNumber(); //numero di step
                //processo che ha creato la traccia 
                if(aStep->GetTrack()->GetCreatorProcess()!=0){                                  //==0 se è una traccia primaria
                    CPName = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
                    CPcode = aStep->GetTrack()->GetCreatorProcess()->GetProcessSubType();   
                    }
                G4ThreeVector vertexpos = aStep->GetTrack()->GetVertexPosition();
                G4ThreeVector vertexmom = aStep->GetTrack()->GetVertexMomentumDirection();
                G4double vertexkinen = aStep->GetTrack()->GetVertexKineticEnergy();
                G4int parentid = aStep->GetTrack()->GetParentID();              //ID della madre
                //G4int parentcode = ->GetDynamicParticle()->GetPDGcode();



                analysisManager->FillNtupleIColumn(0,0,eventNumber);
		        analysisManager->FillNtupleDColumn(0,1,energy);	
                analysisManager->FillNtupleDColumn(0,2,time/ns);
                analysisManager->FillNtupleIColumn(0,3,PDGcode);
                analysisManager->FillNtupleDColumn(0,4,direction.getX());
                analysisManager->FillNtupleDColumn(0,5,direction.getY());
                analysisManager->FillNtupleDColumn(0,6,direction.getZ());
                analysisManager->FillNtupleIColumn(0,7,CPcode);
                analysisManager->FillNtupleDColumn(0,8,vertexpos.getX());
                analysisManager->FillNtupleDColumn(0,9,vertexpos.getY());
                analysisManager->FillNtupleDColumn(0,10,vertexpos.getZ());
                analysisManager->FillNtupleDColumn(0,11,vertexmom.getX());
                analysisManager->FillNtupleDColumn(0,12,vertexmom.getY());
                analysisManager->FillNtupleDColumn(0,13,vertexmom.getZ());
                analysisManager->FillNtupleDColumn(0,14,vertexkinen);
                analysisManager->FillNtupleIColumn(0,15,parentid);
                
                analysisManager->AddNtupleRow(0);
			}
    }


    if(volume->GetName() == "CentralDetector"){
        //cout<<"Got hit"<<endl;
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

            G4double EdepStep = aStep->GetTotalEnergyDeposit();
          //  cout<<"Deposit = "<<EdepStep<<endl;
            fEventAction->addDetectorEdep(0,EdepStep);

            if(aStep->GetPreStepPoint()->GetStepStatus()==fGeomBoundary){       //se entra nel rivelatore
                //fEventAction->addDetectorNh(0,1);
                G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();                 //energia cinetica incidente
                G4ThreeVector direction = aStep->GetPreStepPoint()->GetMomentumDirection();
                G4double time = aStep->GetPreStepPoint()->GetGlobalTime();                      //tempo 
                G4int PDGcode = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();          //PDG code
                //traccia dello step
                G4int CPcode = 0;   
                G4int Nstep = aStep->GetTrack()->GetCurrentStepNumber(); //numero di step
                //processo che ha creato la traccia 
                if(aStep->GetTrack()->GetCreatorProcess()!=0){                                  //==0 se è una traccia primaria
                    CPName = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
                    CPcode = aStep->GetTrack()->GetCreatorProcess()->GetProcessSubType();   
                    }
                G4ThreeVector vertexpos = aStep->GetTrack()->GetVertexPosition();
                G4ThreeVector vertexmom = aStep->GetTrack()->GetVertexMomentumDirection();
                G4double vertexkinen = aStep->GetTrack()->GetVertexKineticEnergy();
                G4int parentid = aStep->GetTrack()->GetParentID();              //ID della madre
                //G4int parentcode = ->GetDynamicParticle()->GetPDGcode();



                analysisManager->FillNtupleIColumn(1,0,eventNumber);
                analysisManager->FillNtupleDColumn(1,1,energy); 
                analysisManager->FillNtupleDColumn(1,2,time/ns);
                analysisManager->FillNtupleIColumn(1,3,PDGcode);
                analysisManager->FillNtupleDColumn(1,4,direction.getX());
                analysisManager->FillNtupleDColumn(1,5,direction.getY());
                analysisManager->FillNtupleDColumn(1,6,direction.getZ());
                analysisManager->FillNtupleIColumn(1,7,CPcode);
                analysisManager->FillNtupleDColumn(1,8,vertexpos.getX());
                analysisManager->FillNtupleDColumn(1,9,vertexpos.getY());
                analysisManager->FillNtupleDColumn(1,10,vertexpos.getZ());
                analysisManager->FillNtupleDColumn(1,11,vertexmom.getX());
                analysisManager->FillNtupleDColumn(1,12,vertexmom.getY());
                analysisManager->FillNtupleDColumn(1,13,vertexmom.getZ());
                analysisManager->FillNtupleDColumn(1,14,vertexkinen);
                analysisManager->FillNtupleIColumn(1,15,parentid);
                
                analysisManager->AddNtupleRow(1);
            }
    }

}





