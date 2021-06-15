#include "SemileptonicTop.h"

using namespace SemileptonicTop;


//************** MC Particle Collection ***************                                                                                                                            
//-----------------------------------------------------

ROOT::VecOps::RVec<int> SemileptonicTop::getMC_daughter(int daughterindex, ROOT::VecOps::RVec<edm4hep::MCParticleData> in,  ROOT::VecOps::RVec<int> ind){
  ROOT::VecOps::RVec<int> result;
  for (size_t i = 0; i < in.size(); ++i) {
    if (daughterindex+1>in.at(i).daughters_end-in.at(i).daughters_begin) {
      result.push_back(-999);
    }
    else {
      result.push_back(ind.at(in.at(i).daughters_begin+daughterindex));
    }
  }
  return result;
}

ROOT::VecOps::RVec<int> SemileptonicTop::getMC_parent(int parentindex, ROOT::VecOps::RVec<edm4hep::MCParticleData> in,  ROOT::VecOps::RVec<int> ind){
  ROOT::VecOps::RVec<int> result;
  for (size_t i = 0; i < in.size(); ++i) {
    if (parentindex+1>in.at(i).parents_end-in.at(i).parents_begin) {
      result.push_back(-999);
    }
    else {
      result.push_back(ind.at(in.at(i).parents_begin+parentindex));
    }
  }
  return result;
}

std::unordered_set<int> SemileptonicTop::selector_finalStates(ROOT::VecOps::RVec<edm4hep::MCParticleData> in, std::unordered_set<int> leptonidx, ROOT::VecOps::RVec<int> RP2MCidx) {
  std::unordered_set<int> result;
  ROOT::VecOps::RVec<int> RP2MCindex{};
  //std::cout << in.size() << " " << leptonidx.size() << " " << RP2MCidx.size() << std::endl;
  for (auto & idx: leptonidx) {
    //std::cout << idx << std::endl;
    RP2MCindex.push_back(RP2MCidx[idx]);
    //std::cout << "RP index for highest energy lepton" << idx << " RP2MC index for highest energy lepton" << RP2MCidx[idx] << std::endl;
  }

  if (RP2MCindex.size() > 1) std::cout << "More than one highest energy lepton found"<< std::endl;
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if (p.generatorStatus!=1) continue;
    if (abs(p.PDG)==12 || abs(p.PDG)==14 || abs(p.PDG)==16) continue;
    if (RP2MCindex.size()>0 && i==RP2MCindex[0]) {
      //std::cout << "Highest energy lepton removed" << std::endl;
      continue;
    }
    result.insert(i);
  }
  return result;
}
/*
std::unordered_set<int> SemileptonicTop::selector_finalStates(ROOT::VecOps::RVec<edm4hep::MCParticleData> in, std::unordered_set<int> leptonidx, ROOT::VecOps::RVec<int> RP2MCidx) {
  std::unordered_set<int> result;

  ROOT::VecOps::RVec<int> RP2MCindex;
  for (auto & idx: leptonidx) {
    RP2MCindex.push_back(RP2MCidx[idx]);
    //std::cout << "RP index for highest energy lepton" << idx << "RP2MC index for highest energy lepton" << RP2MCidx[idx] << std::endl;                                             
  }

  if (RP2MCindex.size() > 1) std::cout << "More than one highest energy lepton found"<< std::endl;

  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if (p.generatorStatus!=1) continue;
    if (abs(p.PDG)==12 || abs(p.PDG)==14 || abs(p.PDG)==16) continue;
    if (i==RP2MCindex[0]) {
      //std::cout << "Highest energy lepton removed" << std::endl;                                                                                                                   
      continue;
    }
    result.insert(i);
  }
  return result;
}
*/

std::unordered_set<int> SemileptonicTop::selector_partons(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  std::unordered_set<int> result;

  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if (p.generatorStatus>80 || p.generatorStatus<70) continue;
    //std::cout << p.generatorStatus << " " << p.PDG << std::endl;                                                                                                                   
    result.insert(i);
  }
  return result;
}


ROOT::VecOps::RVec<edm4hep::MCParticleData> SemileptonicTop::MCParticleSetCreator(ROOT::VecOps::RVec<edm4hep::MCParticleData> in, std::unordered_set<int> idx) {
  ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    std::unordered_set<int>::const_iterator got = idx.find (i);
    if ( got == idx.end() ) continue;
    else result.emplace_back(p);
  }
  return result;
}

ROOT::VecOps::RVec<int> SemileptonicTop::MCParticleSetAssociation(ROOT::VecOps::RVec<edm4hep::MCParticleData> in,std::unordered_set<int> idx) {
  ROOT::VecOps::RVec<int> result{};
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    std::unordered_set<int>::const_iterator got = idx.find (i);
    if ( got == idx.end() ) continue;
    else result.push_back(i);
  }
  return result;
}


//***************** RP Collection *********************                                                                                                                            
//-----------------------------------------------------

std::unordered_set<int> SemileptonicTop::selector_HighestEnergyLepton(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<float> RP2MC_pdg) {
  std::unordered_set<int> result;
  UInt_t nLepton=0;
  float max=0;
  int idx;
  for (size_t i = 0; i < in.size(); ++i) {
    if (abs(RP2MC_pdg[i])==11 || abs(RP2MC_pdg[i])==13) {
      nLepton++;
      auto & p =in[i];
      ROOT::Math::PxPyPzMVector lv(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
      if (nLepton==1) {
        max=lv.E();
        idx=i;
      }
      if (lv.E() > max) {
        max=lv.E();
        idx=i;
      }
    }
  }
  result.insert(idx);
  return result;
}

std::unordered_set<int> SemileptonicTop::selector_rest(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, std::unordered_set<int> idx) {
  std::unordered_set<int> result;
  for (size_t i = 0; i < in.size(); ++i) {
    std::unordered_set<int>::const_iterator got = idx.find (i);
    if ( got == idx.end() ) result.insert(i);
  }
  return result;
}

std::unordered_set<int> SemileptonicTop::vector2selector(std::vector<int> v) {
  std::unordered_set<int> result(v.begin(),v.end());
  return result;
}  

std::vector<edm4hep::ReconstructedParticleData> SemileptonicTop::RPParticleSetCreator(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, std::unordered_set<int> idx) {
  std::vector<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    std::unordered_set<int>::const_iterator got = idx.find (i);
    if ( got == idx.end() ) continue;
    else result.emplace_back(p);
  }
  return result;
}

ROOT::VecOps::RVec<int> SemileptonicTop::RPParticleSetAssociation(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,std::unordered_set<int> idx) {
  ROOT::VecOps::RVec<int> result{};
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    std::unordered_set<int>::const_iterator got = idx.find (i);
    if ( got == idx.end() ) continue;
    else result.push_back(i);
  }
  return result;
}

float SemileptonicTop::RPsetInvariantMass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  float E=0;
  float px=0;
  float py=0;
  float pz=0;
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    ROOT::Math::PxPyPzMVector lv(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    E+=lv.E();
    px+=p.momentum.x;
    py+=p.momentum.y;
    pz+=p.momentum.z;
  }
  float result = sqrt(pow(E,2)-pow(px,2)-pow(py,2)-pow(pz,2));
  return result;

}

float SemileptonicTop::InvMass(ROOT::VecOps::RVec<float> in_px, ROOT::VecOps::RVec<float> in_py, ROOT::VecOps::RVec<float> in_pz, ROOT::VecOps::RVec<float> in_e) {
  float E=0;
  float px=0;
  float py=0;
  float pz=0;
  for (size_t i = 0; i < in_px.size(); ++i) {
    E+=in_e.at(i);
    px+=in_px.at(i);
    py+=in_py.at(i);
    pz+=in_pz.at(i);
  }
  float result = sqrt(pow(E,2)-pow(px,2)-pow(py,2)-pow(pz,2));
  return result;
}

ROOT::VecOps::RVec<float> SemileptonicTop::VertexSignificance(std::vector<std::vector<int>> in, std::vector<edm4hep::ReconstructedParticleData> RPin, ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
  ROOT::VecOps::RVec<float> result;
  for (auto v : in){
    std::unordered_set<int> set_constituents(v.begin(),v.end());
    std::vector<edm4hep::ReconstructedParticleData> RP_constituents = SemileptonicTop::RPParticleSetCreator(RPin, set_constituents);
    VertexingUtils::FCCAnalysesVertex VertexObject = VertexFitterSimple::VertexFitter(1, RP_constituents, tracks); //tracks = EFlowTracks_1
    edm4hep::VertexData Vertex = VertexingUtils::get_VertexData( VertexObject ); //primary vertex, in mm
    float x=Vertex.position.x;
    float y=Vertex.position.y;
    float z=Vertex.position.z;
    float distance2=x*x+y*y+z*z; //distance to IP (0,0,0)  
    float cov_xx = Vertex.covMatrix[0];
    float cov_yx = Vertex.covMatrix[1];
    float cov_yy = Vertex.covMatrix[2];
    float cov_zx = Vertex.covMatrix[3];
    float cov_zy = Vertex.covMatrix[4];
    float cov_zz = Vertex.covMatrix[5];
    float significance=distance2/sqrt(x*x*cov_xx+y*y*cov_yy+z*z*cov_zz+2*(x*y*cov_yx+x*z*cov_zx+y*z*cov_zy));
    result.push_back(significance);
  }
  return result;
}


std::vector<float> SemileptonicTop::alg_sphericity(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
  // ------If number of particles is less than three: special treatment ??? NOT DONE YET
  //-- Compute momentum tensor
  TMatrixD TT(3,3);
  TT.Zero();
  /*
  //--------Test vectors
  std::vector< std::array<float,3> > vect;
  //first vector
  vect.push_back({1, 0, 0});
  //second vector
  vect.push_back({0, 1, 0});
  //third vector
  vect.push_back({0, 0, 1});
  for (std::array<float, 3> p_i : vect) {
    for (int a=0;a<3;a++) {
    for (int b=0;b<3;b++) {
  TT[a][b]+=p_i[a]*p_i[b]; 
}
}
}
  // -------
  */

    for (unsigned int i=0; i<in.size() ; i++) {
      auto & p = in[i];
      float P[3];
      P[0]=p.momentum.x;
      P[1]=p.momentum.y;
      P[2]=p.momentum.z;
      for (int I=0;I<3;I++) {
	for (int K=0;K<I;K++){
	  TT[I][K]+=P[I]*P[K];
	}
      }
    }

    TT[0][1]=TT[1][0];
    TT[0][2]=TT[2][0];
    TT[1][2]=TT[2][1];

    //CALL VSCALE(A,ALPHA,X,N)   X(I)=A(I)*ALPHA   (I=1,..,N)
    float Alpha=1./(TT[0][0]+TT[1][1]+TT[2][2]);
    for (int I=0;I<3;I++) {
      for (int K=0;K<3;K++){
	TT[I][K]*=Alpha;
      }
    }
    const TMatrixDEigen eigen(TT);
    TMatrixD eigenVal = eigen.GetEigenValues();
    TMatrixD eigenVec = eigen.GetEigenVectors();

    std::vector<float> EVAL;
    for (unsigned int i=0; i<3 ; i++) EVAL.push_back(eigenVal[i][i]);

    // Sort the vector in descending order
    std::vector<int> IND(3); //Saving index from sorting eigenvalues
    size_t n(0);
    generate(begin(IND), end(IND), [&]{ return n++; });
    sort(begin(IND),end(IND),[&](int i1, int i2) { return EVAL[i1] > EVAL[i2]; } );

    ///--- DON<C2><B4>T KNOW HOW TO CROSS CHECK FOR IMAGINARY NUMBERS????
    // if (imaginary) continue;
    std::vector<float> major; //major axis (sphericity axis)
    std::vector<float> semimajor; // semi-major axis
    std::vector<float> minor; //minor axis
    for (unsigned int i=0; i<3 ; i++){
      major.push_back(eigenVec[i][IND[0]]);
      semimajor.push_back(eigenVec[i][IND[1]]);
      minor.push_back(eigenVec[i][IND[2]]);
    }
    float sphericity = 1.5*(1-EVAL[IND[0]]);
    float aplanarity = 1.5*EVAL[IND[2]];
    float planarity = EVAL[IND[2]]/EVAL[IND[1]];
    /*    
    std::cout <<"semi-major axis = (" << semimajor[0] << ", " << semimajor[1] << ", " << semimajor[2] << ")" << std::endl;
    std::cout <<"minor axis = (" << minor[0] << ", " << minor[1] << ", " << minor[2] << ")" << std::endl;
    std::cout << "sphericity = " << sphericity << std::endl;
    std::cout << "aplanarity = " << aplanarity << std::endl;
    std::cout << "planarity = " << planarity << std::endl;
  */
    return major;
}
