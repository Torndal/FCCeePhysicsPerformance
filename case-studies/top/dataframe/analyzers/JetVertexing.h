#ifndef  JETVERTEXING_ANALYZERS_H
#define  JETVERTEXING_ANALYZERS_H

#include <cmath>
#include <vector>
#include <unordered_set>

#include "Math/Vector4D.h"
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
#include "SemileptonicTop.h"

/** Jet Vertexing interface.                                                                                                                                                                                                                                                 
This represents a set functions and utilities to loop over number of jets in order to return properties of a vertex.                                                                                                                                                        
*/

namespace JetVertexing{

  struct JetVertex{
    ROOT::VecOps::RVec<float> significance;
    ROOT::VecOps::RVec<float> distance;
    ROOT::VecOps::RVec<float> sigma_d;
    ROOT::VecOps::RVec<float> chi2;
  };

  ROOT::VecOps::RVec<float> get_significance(JetVertex jetvertex);

  ROOT::VecOps::RVec<float> get_distance(JetVertex jetvertex);

  ROOT::VecOps::RVec<float> get_sigma_d(JetVertex jetvertex);

  ROOT::VecOps::RVec<float> get_chi2(JetVertex jetvertex);

  
  //Internal methods
  JetVertex initialise_JetVertex();

  JetVertex build_JetVertex(std::vector<std::vector<int>> in, std::vector<edm4hep::ReconstructedParticleData> RPin, ROOT::VecOps::RVec<edm4hep::TrackState> tracks);
  
}
#endif

