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
/// \file hadronic/Hadr03/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// $Id: PhysicsList.cc 70319 2013-05-29 07:53:09Z gcosmo $

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4IonPhysics.hh"
#include "G4IonINCLXXPhysics.hh"
#include "GammaPhysics.hh"
#include <G4PhysListFactory.hh>
//#include <G4EmLivermorePhysics.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4EmExtraPhysics.hh>
//#include <G4EmStandardPhysics.hh>
#include <G4HadronPhysicsQGSP_BIC_HP.hh>
#include <G4StoppingPhysics.hh>
#include <G4DecayPhysics.hh>
// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Cerenkov.hh"
#include "OpPhysics.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
    :G4VModularPhysicsList()
{
    G4int verb = 1;
    SetVerboseLevel(verb);

    //add new units for cross sections
    //
    new G4UnitDefinition( "mm2/g",  "mm2/g", "Surface/Mass", mm2/g);
    new G4UnitDefinition( "um2/mg", "um2/mg","Surface/Mass", um*um/mg);


    //hadron
    //RegisterPhysics(new G4HadronPhysicsQGSP_BIC_HP(verb));

    // Hadron Elastic scattering
   // RegisterPhysics( new G4HadronElasticPhysics(verb) );

    // Hadron Inelastic Physics
    ////RegisterPhysics( new G4HadronInelasticQBBC(verb));
   // RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(verb));

    ////RegisterPhysics( new G4HadronPhysicsINCLXX(verb));

    // Ion Physics
   // RegisterPhysics( new G4IonPhysics(verb));
    ////RegisterPhysics( new G4IonINCLXXPhysics(ver));
    
    //EM
    //RegisterPhysics(new G4EmLivermorePhysics(verb));
    //RegisterPhysics( new G4EmStandardPhysics(verb) );

    RegisterPhysics( new G4DecayPhysics(verb)); 
    RegisterPhysics( new G4EmStandardPhysics_option4(verb)); 
    RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP(verb)); 
    RegisterPhysics( new G4HadronElasticPhysicsHP(verb)); 
    RegisterPhysics( new G4EmExtraPhysics(verb)); 
    RegisterPhysics( new G4IonPhysics(verb)); 
    RegisterPhysics( new G4StoppingPhysics(verb)); 
    
    // Gamma Physics
    //RegisterPhysics( new GammaPhysics("gamma"));
    RegisterPhysics( new OpPhysics("Cerenkov"));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void PhysicsList::ConstructParticle()
{
    G4BosonConstructor  pBosonConstructor;
    pBosonConstructor.ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();

    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();

    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void PhysicsList::SetCuts()
{
    defaultCutValue = 0.7*mm;
    SetCutValue(0*mm, "proton");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
