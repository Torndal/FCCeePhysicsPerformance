#include "JetVertexing.h"

using namespace JetVertexing;

ROOT::VecOps::RVec<float> JetVertexing::get_significance(JetVertex jetvertex){
  return  jetvertex.significance;
}

ROOT::VecOps::RVec<float> JetVertexing::get_distance(JetVertex jetvertex){
  return  jetvertex.distance;
}

ROOT::VecOps::RVec<float> JetVertexing::get_sigma_d(JetVertex jetvertex){
  return  jetvertex.sigma_d;
}

ROOT::VecOps::RVec<float> JetVertexing::get_chi2(JetVertex jetvertex){
  return  jetvertex.chi2;
}
JetVertex JetVertexing::initialise_JetVertex() {
  JetVertexing::JetVertex result;
  ROOT::VecOps::RVec<float> significance;
  ROOT::VecOps::RVec<float> distance;
  ROOT::VecOps::RVec<float> sigma_d;
  ROOT::VecOps::RVec<float> chi2;

  result.significance = significance;
  result.distance = distance;
  result.sigma_d = sigma_d;
  result.chi2 = chi2;

  return result;
}

JetVertex JetVertexing::build_JetVertex(std::vector<std::vector<int>> in, std::vector<edm4hep::ReconstructedParticleData> RPin, ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  JetVertexing::JetVertex result = JetVertexing::initialise_JetVertex();
  
  for (auto v : in){
    if (v.size() < 2) {
      result.significance.push_back(-999);
      result.distance.push_back(-999);
      result.sigma_d.push_back(-999);
      result.chi2.push_back(-999);
      continue;
    }

    std::unordered_set<int> set_constituents(v.begin(),v.end());
    std::vector<edm4hep::ReconstructedParticleData> RP_constituents = SemileptonicTop::RPParticleSetCreator(RPin, set_constituents);
    VertexingUtils::FCCAnalysesVertex VertexObject = VertexFitterSimple::VertexFitter(1, RP_constituents, tracks); //tracks = EFlowTracks_1                                                                                                                                     
    edm4hep::VertexData Vertex = VertexingUtils::get_VertexData( VertexObject ); //primary vertex, in mm                                                                                                                                                                        
    float x=Vertex.position.x;
    float y=Vertex.position.y;
    float z=Vertex.position.z;
    float distance=sqrt(x*x+y*y+z*z); //distance to IP (0,0,0)
    float cov_xx = Vertex.covMatrix[0];
    float cov_yx = Vertex.covMatrix[1];
    float cov_yy = Vertex.covMatrix[2];
    float cov_zx = Vertex.covMatrix[3];
    float cov_zy = Vertex.covMatrix[4];
    float cov_zz = Vertex.covMatrix[5];
    float sigma_d = sqrt(x*x*cov_xx+y*y*cov_yy+z*z*cov_zz+2*(x*y*cov_yx+x*z*cov_zx+y*z*cov_zy))/distance;
    float significance = distance/sigma_d;
    
    result.significance.push_back(significance);
    result.distance.push_back(distance);
    result.sigma_d.push_back(sigma_d);
    result.chi2.push_back(Vertex.chi2);
    
  }
  return result;
}


