#ifndef  SEMILEPTONICTOP_ANALYZERS_H
#define  SEMILEPTONICTOP_ANALYZERS_H

#include <cmath>
#include <vector>
#include <unordered_set>

#include "Math/Vector4D.h"
#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/Vector3f.h"
#include "edm4hep/Vector3d.h"
#include "edm4hep/Vector2i.h"
#include "edm4hep/TrackState.h"
#include "VertexingUtils.h"
#include "VertexFitterSimple.h"

#include "edm4hep/VertexData.h"
#include "edm4hep/Vertex.h"

#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TVectorT.h"


/** Semileptonic Top interface.
This represents a set functions and utilities to access and perform operations specific to an semileptonic top analysis.
*/

namespace SemileptonicTop{

  //************** MC Particle Collection ***************
  //-----------------------------------------------------
  //EDM4HEP provides parent and daughter history for the particles in the MC Collection. This corresponds to the mother and daughter columns in the Pythia output.
  //Note that EDM4HEP also saves the history for ISR particles unlike Pythia
  //Returns the index of the daugther particle in the MCParticle collection for the index given for the column
  ROOT::VecOps::RVec<int> getMC_daughter(int daughterindex, ROOT::VecOps::RVec<edm4hep::MCParticleData> in,  ROOT::VecOps::RVec<int> ind);
  //Returns the index of the parent particle in the MCParticle collection for the index given for the column
  ROOT::VecOps::RVec<int> getMC_parent(int parentindex, ROOT::VecOps::RVec<edm4hep::MCParticleData> in,  ROOT::VecOps::RVec<int> ind);

  //Select final state particles except the lepton corresponding to the highest energy lepton in the RP collection
  std::unordered_set<int> selector_finalStates(ROOT::VecOps::RVec<edm4hep::MCParticleData> in, std::unordered_set<int> leptonidx, ROOT::VecOps::RVec<int> RP2MCidx);

  //Select particles with status codes between 71-79 which corresponds to partons in preparation of hadronization process
  std::unordered_set<int> selector_partons(ROOT::VecOps::RVec<edm4hep::MCParticleData> in);

  //Particle Set Creator using the selectors
  ROOT::VecOps::RVec<edm4hep::MCParticleData> MCParticleSetCreator(ROOT::VecOps::RVec<edm4hep::MCParticleData> in, std::unordered_set<int> idx);

  //Provides associations between the subset of particles and the full MC collection. 
  ROOT::VecOps::RVec<int> MCParticleSetAssociation(ROOT::VecOps::RVec<edm4hep::MCParticleData> in,std::unordered_set<int> idx);


  //***************** RP Collection *********************
  //-----------------------------------------------------
  //return indices for a particle set of highest energy leptons
  std::unordered_set<int> selector_HighestEnergyLepton(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<float> RP2MC_pdg);
  
  //return indices for the remaining particle set
  std::unordered_set<int> selector_rest(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, std::unordered_set<int> idx);

  //convert vector to unordered_set 
  std::unordered_set<int> vector2selector(std::vector<int> v);
  
  //return particle set from some selector
  std::vector<edm4hep::ReconstructedParticleData> RPParticleSetCreator(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, std::unordered_set<int> idx);
  
  //return particle set index in full RP collection as a vector
  ROOT::VecOps::RVec<int> RPParticleSetAssociation(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,std::unordered_set<int> idx);

  //return invariant mass of the input collection
  float RPsetInvariantMass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  ROOT::VecOps::RVec<float> VertexSignificance(std::vector<std::vector<int>> in, std::vector<edm4hep::ReconstructedParticleData> RPin, ROOT::VecOps::RVec<edm4hep::TrackState> tracks);
  
  //return sphericity based on ALEPH code
  std::vector<float> alg_sphericity(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

}
#endif
