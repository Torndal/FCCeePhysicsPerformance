import sys
import ROOT

print ("Load cxx analyzers ... ",)
ROOT.gSystem.Load("libedm4hep")
ROOT.gSystem.Load("libpodio")
ROOT.gSystem.Load("libFCCAnalyses")
ROOT.gSystem.Load("libfastjet")
ROOT.gSystem.Load("libFCCAnalysesTop")
ROOT.gErrorIgnoreLevel = ROOT.kFatal
_edm  = ROOT.edm4hep.ReconstructedParticleData()
_pod  = ROOT.podio.ObjectID()
_fcc  = ROOT.dummyLoader
_top  = ROOT.dummyLoaderTop

print ('edm4hep  ',_edm)
print ('podio    ',_pod)
print ('fccana   ',_fcc)
print ('fcctop   ',_top)

class analysis():

    #__________________________________________________________
    def __init__(self, inputlist, outname, ncpu):
        self.outname = outname
        if ".root" not in outname:
            self.outname+=".root"

        ROOT.ROOT.EnableImplicitMT(ncpu)

        self.df = ROOT.RDataFrame("events", inputlist)
        print (" done")
    #__________________________________________________________
    def run(self):
        #match=ROOT.getRP2MC_p_func()
        string_vec = ROOT.std.vector('string')()
        string_vec.push_back('MCRecoAssociations#0.index')
        string_vec.push_back('MCRecoAssociations#1.index')
        string_vec.push_back('ReconstructedParticles')
        string_vec.push_back('Particle')

        df2 = (self.df
               .Define("MC_px",         "MCParticle::get_px(Particle)")
               .Define("MC_py",         "MCParticle::get_py(Particle)")
               .Define("MC_pz",         "MCParticle::get_pz(Particle)")
               .Define("MC_p",          "MCParticle::get_p(Particle)")
               .Define("MC_pdg",        "MCParticle::get_pdg(Particle)")
               .Define("MC_charge",     "MCParticle::get_charge(Particle)")
               .Define("MC_mass",       "MCParticle::get_mass(Particle)")
               .Define("MC_e",          "MCParticle::get_e(Particle)")
               .Define("MC_status",     "MCParticle::get_genStatus(Particle)")

               .Define("RP_p",          "ReconstructedParticle::get_p(ReconstructedParticles)")
               .Define("RP_px",         "ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("RP_py",         "ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("RP_pz",         "ReconstructedParticle::get_pz(ReconstructedParticles)")
               .Define("RP_charge",     "ReconstructedParticle::get_charge(ReconstructedParticles)")
               .Define("RP_mass",       "ReconstructedParticle::get_mass(ReconstructedParticles)")
               .Define("RP_e",          "ReconstructedParticle::get_e(ReconstructedParticles)")
               
               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
               .Alias("Particle0", "Particle#0.index")
               .Alias("Particle1", "Particle#1.index")

               .Define('RPMC_index',    "ReconstructedParticle2MC::getRP2MC_index(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles)")
               .Define('RPMC_pdg',      "ReconstructedParticle2MC::getRP2MC_pdg(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               .Define('RPMC_parentindex', "MCParticle::get_parentid(RPMC_index,Particle, Particle0)")

               .Define('MC_daughter1', "SemileptonicTop::getMC_daughter(0,Particle,Particle1)")
               .Define('MC_daughter2', "SemileptonicTop::getMC_daughter(1,Particle,Particle1)")
               .Define('MC_parent1', "SemileptonicTop::getMC_parent(0,Particle,Particle0)")
               .Define('MC_parent2', "SemileptonicTop::getMC_parent(1,Particle,Particle0)")

               
               .Define("MET_p",          "ReconstructedParticle::get_p(MissingET)")
               .Define("MET_px",         "ReconstructedParticle::get_px(MissingET)")
               .Define("MET_py",         "ReconstructedParticle::get_py(MissingET)")
               .Define("MET_pz",         "ReconstructedParticle::get_pz(MissingET)")
               .Define("MET_charge",     "ReconstructedParticle::get_charge(MissingET)")
               .Define("MET_mass",       "ReconstructedParticle::get_mass(MissingET)")
               .Define("MET_e",          "ReconstructedParticle::get_e(MissingET)")

               .Define("Set_lepton",     "SemileptonicTop::selector_HighestEnergyLepton(ReconstructedParticles, RPMC_pdg)")
               .Define("RPlepton",       "SemileptonicTop::RPParticleSetCreator(ReconstructedParticles, Set_lepton)")
               .Define("RPlepton_p",     "ReconstructedParticle::get_p(RPlepton)")
               .Define("RPleptonMET",    "ReconstructedParticle::merge(RPlepton,MissingET)")
               .Define("RPleptonMET_invmass","SemileptonicTop::RPsetInvariantMass(RPleptonMET)")
               
               .Define("Set_rest",       "SemileptonicTop::selector_rest(ReconstructedParticles,Set_lepton)")
               .Define("RPrest",         "SemileptonicTop::RPParticleSetCreator(ReconstructedParticles, Set_rest)")
               .Define("RPrest_association", "SemileptonicTop::RPParticleSetAssociation(ReconstructedParticles, Set_rest)")
               .Define("RPrest_invmass", "SemileptonicTop::RPsetInvariantMass(RPrest)")
               .Define("RPrest_px",      "ReconstructedParticle::get_px(RPrest)")
               .Define("RPrest_py",      "ReconstructedParticle::get_py(RPrest)")
               .Define("RPrest_pz",      "ReconstructedParticle::get_pz(RPrest)")
               .Define("RPrest_charge",      "ReconstructedParticle::get_charge(RPrest)")
               .Define("RPrest_e",          "ReconstructedParticle::get_e(RPrest)")
               
               .Define("Set_finalStates", "SemileptonicTop::selector_finalStates(Particle, Set_lepton, RPMC_index)")
               .Define("Particle_finalStates", "SemileptonicTop::MCParticleSetCreator(Particle, Set_finalStates)")
               .Define("MCfinal_association", "SemileptonicTop::MCParticleSetAssociation(Particle, Set_finalStates)")
               .Define("MCfinal_px", "MCParticle::get_px(Particle_finalStates)")
               .Define("MCfinal_py", "MCParticle::get_py(Particle_finalStates)")
               .Define("MCfinal_pz", "MCParticle::get_pz(Particle_finalStates)")
               .Define("MCfinal_e",  "MCParticle::get_e(Particle_finalStates)")

               .Define("Set_partons", "SemileptonicTop::selector_partons(Particle)")
               .Define("Particle_partons", "SemileptonicTop::MCParticleSetCreator(Particle, Set_partons)")
               .Define("MCparton_association", "SemileptonicTop::MCParticleSetAssociation(Particle, Set_partons)")
               .Define("MCparton_px", "MCParticle::get_px(Particle_partons)")
               .Define("MCparton_py", "MCParticle::get_py(Particle_partons)")
               .Define("MCparton_pz", "MCParticle::get_pz(Particle_partons)")
               .Define("MCparton_e",  "MCParticle::get_e(Particle_partons)")


               .Define('EVT_thrust',      'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               .Define('EVT_thrust_val',  'EVT_thrust.at(0)')
               .Define('EVT_thrust_x',    'EVT_thrust.at(1)')
               .Define('EVT_thrust_x_err','EVT_thrust.at(2)')
               .Define('EVT_thrust_y',    'EVT_thrust.at(3)')
               .Define('EVT_thrust_y_err','EVT_thrust.at(4)')
               .Define('EVT_thrust_z',    'EVT_thrust.at(5)')
               .Define('EVT_thrust_z_err','EVT_thrust.at(6)')
               
               .Define('EVT_sphericity',     'Algorithms::minimize_sphericity("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               .Define('EVT_sphericity_x',   "EVT_sphericity.at(0)")
               .Define('EVT_sphericity_y',   "EVT_sphericity.at(1)")
               .Define('EVT_sphericity_z',   "EVT_sphericity.at(2)")
               .Define('EVT_sphericity_val', "EVT_sphericity.at(3)")

               .Define('EVTrest_thrust',     'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               .Define('EVTrest_thrust_val',  'EVTrest_thrust.at(0)')
               .Define('EVTrest_thrust_x',    'EVTrest_thrust.at(1)')
               .Define('EVTrest_thrust_x_err','EVTrest_thrust.at(2)')
               .Define('EVTrest_thrust_y',    'EVTrest_thrust.at(3)')
               .Define('EVTrest_thrust_y_err','EVTrest_thrust.at(4)')
               .Define('EVTrest_thrust_z',    'EVTrest_thrust.at(5)')
               .Define('EVTrest_thrust_z_err','EVTrest_thrust.at(6)')

               .Define("AlgoSpher", "SemileptonicTop::alg_sphericity(ReconstructedParticles)")

               #build pseudo jets from RP set (excluding HE lepton), MC set (exluding RP2MC HE lepton), and parton set
               .Define("reco_jets",    "JetClusteringUtils::set_pseudoJets(RPrest_px, RPrest_py, RPrest_pz, RPrest_e)")
               .Define("particle_jets", "JetClusteringUtils::set_pseudoJets(MCfinal_px, MCfinal_py, MCfinal_pz, MCfinal_e)")
               .Define("parton_jets",  "JetClusteringUtils::set_pseudoJets(MCparton_px, MCparton_py, MCparton_pz, MCparton_e)")


               ################## kT ###################
               #run jet clustering with reconstructed particles. kt_algorithm, R=1.0, exclusive clustering, exactly 4 jets, sorted by E, E-scheme
               .Define("FCCAnalysesRecoJets_kt", "JetClustering::clustering_kt(1.0, 2, 4, 1, 0)(reco_jets)")
               #get the jets out of the struct
               .Define("recojets_kt",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesRecoJets_kt)")
               #get the jets constituents out of the struct
               .Define("recojetconstituents_kt","JetClusteringUtils::get_constituents(FCCAnalysesRecoJets_kt)")
               #get some variables
               .Define("recojets_kt_px",        "JetClusteringUtils::get_px(recojets_kt)")
               .Define("recojets_kt_py",        "JetClusteringUtils::get_py(recojets_kt)")
               .Define("recojets_kt_pz",        "JetClusteringUtils::get_pz(recojets_kt)")
               .Define("recojets_kt_e",         "JetClusteringUtils::get_e(recojets_kt)")

               .Define("FCCAnalysesParticleJets_kt", "JetClustering::clustering_kt(1.0, 2, 4, 1, 0)(particle_jets)")
               .Define("particlejets_kt",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesParticleJets_kt)")
               .Define("particlejetconstituents_kt","JetClusteringUtils::get_constituents(FCCAnalysesParticleJets_kt)")
               .Define("particlejets_kt_px",        "JetClusteringUtils::get_px(particlejets_kt)")
               .Define("particlejets_kt_py",        "JetClusteringUtils::get_py(particlejets_kt)")
               .Define("particlejets_kt_pz",        "JetClusteringUtils::get_pz(particlejets_kt)")
               .Define("particlejets_kt_e",        "JetClusteringUtils::get_e(particlejets_kt)")

               .Define("FCCAnalysesPartonJets_kt", "JetClustering::clustering_kt(1.0, 2, 4, 1, 0)(parton_jets)")
               .Define("partonjets_kt",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesPartonJets_kt)")
               .Define("partonjetconstituents_kt","JetClusteringUtils::get_constituents(FCCAnalysesPartonJets_kt)")
               .Define("partonjets_kt_px",        "JetClusteringUtils::get_px(partonjets_kt)")
               .Define("partonjets_kt_py",        "JetClusteringUtils::get_py(partonjets_kt)")
               .Define("partonjets_kt_pz",        "JetClusteringUtils::get_pz(partonjets_kt)")
               .Define("partonjets_kt_e",        "JetClusteringUtils::get_e(partonjets_kt)")

               ################## Durham ###################
               #run jet clustering with reconstructed particles. Durham_algorithm, exclusive clustering, exactly 4 jets, sorted by E, E-scheme  
               .Define("FCCAnalysesRecoJets_durham", "JetClustering::clustering_ee_kt(2, 4, 1, 0)(reco_jets)")
               .Define("recojets_durham",            "JetClusteringUtils::get_pseudoJets(FCCAnalysesRecoJets_durham)")
               .Define("recojetconstituents_durham", "JetClusteringUtils::get_constituents(FCCAnalysesRecoJets_durham)")
               .Define("recojets_durham_dmin",       "JetClusteringUtils::get_dmerge(FCCAnalysesRecoJets_durham)")
               .Define("recojetymerge_durham",       "JetClusteringUtils::get_ymerge(FCCAnalysesRecoJets_durham)")
               .Define("recojets_durham_px",         "JetClusteringUtils::get_px(recojets_durham)")
               .Define("recojets_durham_py",         "JetClusteringUtils::get_py(recojets_durham)")
               .Define("recojets_durham_pz",         "JetClusteringUtils::get_pz(recojets_durham)")
               .Define("recojets_durham_e",          "JetClusteringUtils::get_e(recojets_durham)")
               .Define("recojets_durham_flavour",    "JetTaggingUtils::get_flavour(recojets_durham, Particle)")
               .Define("recojets_durham_btag",       "JetTaggingUtils::get_btag(recojets_durham_flavour, 1.00)")
               
               .Define("FCCAnalysesParticleJets_durham", "JetClustering::clustering_ee_kt(2, 4, 1, 0)(particle_jets)")
               .Define("particlejets_durham",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesParticleJets_durham)")
               .Define("particlejetconstituents_durham","JetClusteringUtils::get_constituents(FCCAnalysesParticleJets_durham)")
               .Define("particlejets_durham_px",        "JetClusteringUtils::get_px(particlejets_durham)")
               .Define("particlejets_durham_py",        "JetClusteringUtils::get_py(particlejets_durham)")
               .Define("particlejets_durham_pz",        "JetClusteringUtils::get_pz(particlejets_durham)")
               .Define("particlejets_durham_e",        "JetClusteringUtils::get_e(particlejets_durham)")
               .Define("particlejets_durham_flavour",    "JetTaggingUtils::get_flavour(particlejets_durham, Particle)")
               .Define("particlejets_durham_btag",       "JetTaggingUtils::get_btag(particlejets_durham_flavour, 1.00)")


               .Define("FCCAnalysesPartonJets_durham", "JetClustering::clustering_ee_kt(2, 4, 1, 0)(parton_jets)")
               .Define("partonjets_durham",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesPartonJets_durham)")
               .Define("partonjetconstituents_durham","JetClusteringUtils::get_constituents(FCCAnalysesPartonJets_durham)")
               .Define("partonjets_durham_px",        "JetClusteringUtils::get_px(partonjets_durham)")
               .Define("partonjets_durham_py",        "JetClusteringUtils::get_py(partonjets_durham)")
               .Define("partonjets_durham_pz",        "JetClusteringUtils::get_pz(partonjets_durham)")
               .Define("partonjets_durham_e",        "JetClusteringUtils::get_e(partonjets_durham)")
               .Define("partonjets_durham_flavour",    "JetTaggingUtils::get_flavour(partonjets_durham, Particle)")
               .Define("partonjets_durham_btag",       "JetTaggingUtils::get_btag(partonjets_durham_flavour, 1.00)")

               #run jet clustering with reconstructed particles. Durham_algorithm, exclusive clustering, exactly 4 jets, sorted by E, E0-scheme  
               .Define("FCCAnalysesRecoJets_durhamE0", "JetClustering::clustering_ee_kt(2, 4, 1, 10)(reco_jets)")
               .Define("recojets_durhamE0",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesRecoJets_durhamE0)")
               .Define("recojetconstituents_durhamE0","JetClusteringUtils::get_constituents(FCCAnalysesRecoJets_durhamE0)")
               .Define("recojets_durhamE0_px",        "JetClusteringUtils::get_px(recojets_durhamE0)")
               .Define("recojets_durhamE0_py",        "JetClusteringUtils::get_py(recojets_durhamE0)")
               .Define("recojets_durhamE0_pz",        "JetClusteringUtils::get_pz(recojets_durhamE0)")
               .Define("recojets_durhamE0_e",        "JetClusteringUtils::get_e(recojets_durhamE0)")
               .Define("recojets_durhamE0_flavour",    "JetTaggingUtils::get_flavour(recojets_durhamE0, Particle)")
               .Define("recojets_durhamE0_btag",       "JetTaggingUtils::get_btag(recojets_durhamE0_flavour, 1.00)")
               
               .Define("FCCAnalysesParticleJets_durhamE0", "JetClustering::clustering_ee_kt(2, 4, 1, 10)(particle_jets)")
               .Define("particlejets_durhamE0",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesParticleJets_durhamE0)")
               .Define("particlejetconstituents_durhamE0","JetClusteringUtils::get_constituents(FCCAnalysesParticleJets_durhamE0)")
               .Define("particlejets_durhamE0_px",        "JetClusteringUtils::get_px(particlejets_durhamE0)")
               .Define("particlejets_durhamE0_py",        "JetClusteringUtils::get_py(particlejets_durhamE0)")
               .Define("particlejets_durhamE0_pz",        "JetClusteringUtils::get_pz(particlejets_durhamE0)")
               .Define("particlejets_durhamE0_e",        "JetClusteringUtils::get_e(particlejets_durhamE0)")

               #run jet clustering with reconstructed particles. Durham_algorithm, exclusive clustering, exactly 4 jets, sorted by E, p-scheme  
               .Define("FCCAnalysesRecoJets_durhamp", "JetClustering::clustering_ee_kt(2, 4, 1,11)(reco_jets)")
               .Define("recojets_durhamp",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesRecoJets_durhamp)")
               .Define("recojetconstituents_durhamp","JetClusteringUtils::get_constituents(FCCAnalysesRecoJets_durhamp)")
               .Define("recojets_durhamp_px",        "JetClusteringUtils::get_px(recojets_durhamp)")
               .Define("recojets_durhamp_py",        "JetClusteringUtils::get_py(recojets_durhamp)")
               .Define("recojets_durhamp_pz",        "JetClusteringUtils::get_pz(recojets_durhamp)")
               .Define("recojets_durhamp_e",        "JetClusteringUtils::get_e(recojets_durhamp)")
               .Define("recojets_durhamp_flavour",    "JetTaggingUtils::get_flavour(recojets_durhamp, Particle)")
               .Define("recojets_durhamp_btag",       "JetTaggingUtils::get_btag(recojets_durhamp_flavour, 1.00)")
               
               .Define("FCCAnalysesParticleJets_durhamp", "JetClustering::clustering_ee_kt(2, 4, 1, 11)(particle_jets)")
               .Define("particlejets_durhamp",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesParticleJets_durhamp)")
               .Define("particlejetconstituents_durhamp","JetClusteringUtils::get_constituents(FCCAnalysesParticleJets_durhamp)")
               .Define("particlejets_durhamp_px",        "JetClusteringUtils::get_px(particlejets_durhamp)")
               .Define("particlejets_durhamp_py",        "JetClusteringUtils::get_py(particlejets_durhamp)")
               .Define("particlejets_durhamp_pz",        "JetClusteringUtils::get_pz(particlejets_durhamp)")
               .Define("particlejets_durhamp_e",        "JetClusteringUtils::get_e(particlejets_durhamp)")

               #------------ Vertexing ---------------------
               #.Define("recojetconstituents_durham0", "recojetconstituents_durham.at(0)")
               #.Define("Set_constituents", "SemileptonicTop::vector2selector(recojetconstituents_durham0)")
               #constituents point to the subset of particle that's needed as input to get the right association
               #.Define("RP_constituents",         "SemileptonicTop::RPParticleSetCreator(RPrest, Set_constituents)")
               #.Define("VertexObject",  "VertexFitterSimple::VertexFitter( 1, RP_constituents, EFlowTrack_1 )")
               #.Define("Vertex",        "VertexingUtils::get_VertexData( VertexObject )")    # primary vertex, in mm
               #.Define("Vertex_x", "Vertex.position.x")
               #.Define("Vertex_y", "Vertex.position.y")
               #.Define("Vertex_z", "Vertex.position.z")
               #.Define("covMatrix_xx", "Vertex.covMatrix[0]")
               #.Define("covMatrix_yx", "Vertex.covMatrix[1]")
               #.Define("covMatrix_zx", "Vertex.covMatrix[2]")
               #.Define("covMatrix_yy", "Vertex.covMatrix[3]")
               #.Define("covMatrix_zy", "Vertex.covMatrix[4]")
               #.Define("covMatrix_zz", "Vertex.covMatrix[5]")

               .Define("Significance_reco", "SemileptonicTop::VertexSignificance(recojetconstituents_durham, RPrest, EFlowTrack_1)")
               .Define("Significance_particle", "SemileptonicTop::VertexSignificance(particlejetconstituents_durham, RPrest, EFlowTrack_1)")
               .Define("Significance_parton", "SemileptonicTop::VertexSignificance(partonjetconstituents_durham, RPrest, EFlowTrack_1)")
               
               ################## ee anti-kT ###################
               #run jet clustering with reconstructed particles. ee_anti-kT_algorithm, exclusive clustering, exactly 4 jets, sorted by E, E-scheme 
               #Even though in FastJet User Manual they advise against the use of exclusive jets in the context of the anti-kT algorithm,
               #because of the lack of physically meaningful hierarchy in the clustering sequence 
               .Define("FCCAnalysesRecoJets_ee_antikt", "JetClustering::clustering_ee_genkt(1.0, 2, 4, 1, 0, -1)(reco_jets)")
               .Define("recojets_ee_antikt",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesRecoJets_ee_antikt)")
               .Define("recojetconstituents_ee_antikt","JetClusteringUtils::get_constituents(FCCAnalysesRecoJets_ee_antikt)")
               .Define("recojets_ee_antikt_px",        "JetClusteringUtils::get_px(recojets_ee_antikt)")
               .Define("recojets_ee_antikt_py",        "JetClusteringUtils::get_py(recojets_ee_antikt)")
               .Define("recojets_ee_antikt_pz",        "JetClusteringUtils::get_pz(recojets_ee_antikt)")
               .Define("recojets_ee_antikt_e",        "JetClusteringUtils::get_e(recojets_ee_antikt)")

               .Define("FCCAnalysesParticleJets_ee_antikt", "JetClustering::clustering_ee_genkt(1.0, 2, 4, 1, 0, -1)(particle_jets)")
               .Define("particlejets_ee_antikt",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesParticleJets_ee_antikt)")
               .Define("particlejetconstituents_ee_antikt","JetClusteringUtils::get_constituents(FCCAnalysesParticleJets_ee_antikt)")
               .Define("particlejets_ee_antikt_px",        "JetClusteringUtils::get_px(particlejets_ee_antikt)")
               .Define("particlejets_ee_antikt_py",        "JetClusteringUtils::get_py(particlejets_ee_antikt)")
               .Define("particlejets_ee_antikt_pz",        "JetClusteringUtils::get_pz(particlejets_ee_antikt)")
               .Define("particlejets_ee_antikt_e",        "JetClusteringUtils::get_e(particlejets_ee_antikt)")

               .Define("FCCAnalysesPartonJets_ee_antikt", "JetClustering::clustering_ee_genkt(1.0, 2, 4, 1, 0, -1)(parton_jets)")
               .Define("partonjets_ee_antikt",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesPartonJets_ee_antikt)")
               .Define("partonjetconstituents_ee_antikt","JetClusteringUtils::get_constituents(FCCAnalysesPartonJets_ee_antikt)")
               .Define("partonjets_ee_antikt_px",        "JetClusteringUtils::get_px(partonjets_ee_antikt)")
               .Define("partonjets_ee_antikt_py",        "JetClusteringUtils::get_py(partonjets_ee_antikt)")
               .Define("partonjets_ee_antikt_pz",        "JetClusteringUtils::get_pz(partonjets_ee_antikt)")
               .Define("partonjets_ee_antikt_e",        "JetClusteringUtils::get_e(partonjets_ee_antikt)")

               ################## Cambridge ###################
               #run jet clustering with reconstructed particles. ee_cambridge_algorithm, exclusive clustering, exactly 4 jets, sorted by E, E-scheme 
               .Define("FCCAnalysesRecoJets_cambridge", "JetClustering::clustering_ee_genkt(1.0, 2, 4, 1, 0, 0)(reco_jets)")
               .Define("recojets_cambridge",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesRecoJets_cambridge)")
               .Define("recojetconstituents_cambridge","JetClusteringUtils::get_constituents(FCCAnalysesRecoJets_cambridge)")
               .Define("recojets_cambridge_px",        "JetClusteringUtils::get_px(recojets_cambridge)")
               .Define("recojets_cambridge_py",        "JetClusteringUtils::get_py(recojets_cambridge)")
               .Define("recojets_cambridge_pz",        "JetClusteringUtils::get_pz(recojets_cambridge)")
               .Define("recojets_cambridge_e",        "JetClusteringUtils::get_e(recojets_cambridge)")

               .Define("FCCAnalysesParticleJets_cambridge", "JetClustering::clustering_ee_genkt(1.0, 2, 4, 1, 0, 0)(particle_jets)")
               .Define("particlejets_cambridge",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesParticleJets_cambridge)")
               .Define("particlejetconstituents_cambridge","JetClusteringUtils::get_constituents(FCCAnalysesParticleJets_cambridge)")
               .Define("particlejets_cambridge_px",        "JetClusteringUtils::get_px(particlejets_cambridge)")
               .Define("particlejets_cambridge_py",        "JetClusteringUtils::get_py(particlejets_cambridge)")
               .Define("particlejets_cambridge_pz",        "JetClusteringUtils::get_pz(particlejets_cambridge)")
               .Define("particlejets_cambridge_e",        "JetClusteringUtils::get_e(particlejets_cambridge)")

               .Define("FCCAnalysesPartonJets_cambridge", "JetClustering::clustering_ee_genkt(1.0, 2, 4, 1, 0, 0)(parton_jets)")
               .Define("partonjets_cambridge",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesPartonJets_cambridge)")
               .Define("partonjetconstituents_cambridge","JetClusteringUtils::get_constituents(FCCAnalysesPartonJets_cambridge)")
               .Define("partonjets_cambridge_px",        "JetClusteringUtils::get_px(partonjets_cambridge)")
               .Define("partonjets_cambridge_py",        "JetClusteringUtils::get_py(partonjets_cambridge)")
               .Define("partonjets_cambridge_pz",        "JetClusteringUtils::get_pz(partonjets_cambridge)")
               .Define("partonjets_cambridge_e",        "JetClusteringUtils::get_e(partonjets_cambridge)")
               
               ################## Jade ###################
               #run jet clustering with reconstructed particles. Jade_algorithm, exclusive clustering, exactly 4 jets, sorted by E, E-scheme 
               .Define("FCCAnalysesRecoJets_jade", "JetClustering::clustering_jade(1.0, 2, 4, 1, 0)(reco_jets)")
               .Define("recojets_jade",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesRecoJets_jade)")
               .Define("recojetconstituents_jade","JetClusteringUtils::get_constituents(FCCAnalysesRecoJets_jade)")
               .Define("recojets_jade_px",        "JetClusteringUtils::get_px(recojets_jade)")
               .Define("recojets_jade_py",        "JetClusteringUtils::get_py(recojets_jade)")
               .Define("recojets_jade_pz",        "JetClusteringUtils::get_pz(recojets_jade)")
               .Define("recojets_jade_e",        "JetClusteringUtils::get_e(recojets_jade)")

               .Define("FCCAnalysesParticleJets_jade", "JetClustering::clustering_jade(1.0, 2, 4, 1, 0)(particle_jets)")
               .Define("particlejets_jade",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesParticleJets_jade)")
               .Define("particlejetconstituents_jade","JetClusteringUtils::get_constituents(FCCAnalysesParticleJets_jade)")
               .Define("particlejets_jade_px",        "JetClusteringUtils::get_px(particlejets_jade)")
               .Define("particlejets_jade_py",        "JetClusteringUtils::get_py(particlejets_jade)")
               .Define("particlejets_jade_pz",        "JetClusteringUtils::get_pz(particlejets_jade)")
               .Define("particlejets_jade_e",        "JetClusteringUtils::get_e(particlejets_jade)")

               .Define("FCCAnalysesPartonJets_jade", "JetClustering::clustering_jade(1.0, 2, 4, 1, 0)(parton_jets)")
               .Define("partonjets_jade",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesPartonJets_jade)")
               .Define("partonjetconstituents_jade","JetClusteringUtils::get_constituents(FCCAnalysesPartonJets_jade)")
               .Define("partonjets_jade_px",        "JetClusteringUtils::get_px(partonjets_jade)")
               .Define("partonjets_jade_py",        "JetClusteringUtils::get_py(partonjets_jade)")
               .Define("partonjets_jade_pz",        "JetClusteringUtils::get_pz(partonjets_jade)")
               .Define("partonjets_jade_e",        "JetClusteringUtils::get_e(partonjets_jade)")

               ################## Valencia ###################
               #run jet clustering with reconstructed particles. Valencia_algorithm, exclusive clustering, exactly 4 jets, sorted by E, E-scheme 
               .Define("FCCAnalysesRecoJets_valencia", "JetClustering::clustering_valencia(1.0, 2, 4, 1, 0, 1., 1.)(reco_jets)")
               .Define("recojets_valencia",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesRecoJets_valencia)")
               .Define("recojetconstituents_valencia","JetClusteringUtils::get_constituents(FCCAnalysesRecoJets_valencia)")
               .Define("recojets_valencia_px",        "JetClusteringUtils::get_px(recojets_valencia)")
               .Define("recojets_valencia_py",        "JetClusteringUtils::get_py(recojets_valencia)")
               .Define("recojets_valencia_pz",        "JetClusteringUtils::get_pz(recojets_valencia)")
               .Define("recojets_valencia_e",        "JetClusteringUtils::get_e(recojets_valencia)")

               .Define("FCCAnalysesParticleJets_valencia", "JetClustering::clustering_valencia(1.0, 2, 4, 1, 0, 1., 1.)(particle_jets)")
               .Define("particlejets_valencia",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesParticleJets_valencia)")
               .Define("particlejetconstituents_valencia","JetClusteringUtils::get_constituents(FCCAnalysesParticleJets_valencia)")
               .Define("particlejets_valencia_px",        "JetClusteringUtils::get_px(particlejets_valencia)")
               .Define("particlejets_valencia_py",        "JetClusteringUtils::get_py(particlejets_valencia)")
               .Define("particlejets_valencia_pz",        "JetClusteringUtils::get_pz(particlejets_valencia)")
               .Define("particlejets_valencia_e",        "JetClusteringUtils::get_e(particlejets_valencia)")

               .Define("FCCAnalysesPartonJets_valencia", "JetClustering::clustering_valencia(1.0, 2, 4, 1, 0, 1., 1.)(parton_jets)")
               .Define("partonjets_valencia",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesPartonJets_valencia)")
               .Define("partonjetconstituents_valencia","JetClusteringUtils::get_constituents(FCCAnalysesPartonJets_valencia)")
               .Define("partonjets_valencia_px",        "JetClusteringUtils::get_px(partonjets_valencia)")
               .Define("partonjets_valencia_py",        "JetClusteringUtils::get_py(partonjets_valencia)")
               .Define("partonjets_valencia_pz",        "JetClusteringUtils::get_pz(partonjets_valencia)")
               .Define("partonjets_valencia_e",        "JetClusteringUtils::get_e(partonjets_valencia)")

               )

        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [
                "MC_px",
                "MC_py",
                "MC_pz",
                "MC_p",
                "MC_pdg",
                "MC_charge",
                "MC_mass",
                "MC_e",
                "MC_status",
               
                "RP_p",
                "RP_px",
                "RP_py",
                "RP_pz",
                "RP_charge",
                "RP_mass",
                "RP_e",

                "RPMC_pdg",
                "RPMC_index",
                "RPMC_parentindex",

                "MC_daughter1",
                "MC_daughter2",
                "MC_parent1",
                "MC_parent2",
                
                "MET_p",
                "MET_px",
                "MET_py",
                "MET_pz",
                "MET_charge",
                "MET_mass",
                "MET_e",

                "RPlepton_p",
                "RPleptonMET_invmass",
                "RPrest_invmass",
                "RPrest_association",
                "MCfinal_association",
                "MCparton_association",
                
                "EVT_thrust_val",
                "EVTrest_thrust_val",
                
                #"AlgoSpher",

                "recojets_kt_px",
                "recojets_kt_py",
                "recojets_kt_pz",
                "recojets_kt_e",
                "recojetconstituents_kt",
                "particlejets_kt_px",
                "particlejets_kt_py",
                "particlejets_kt_pz",
                "particlejets_kt_e",
                "particlejetconstituents_kt",
                "partonjets_kt_px",
                "partonjets_kt_py",
                "partonjets_kt_pz",
                "partonjets_kt_e",
                "partonjetconstituents_kt",

                "recojets_durham_px",
                "recojets_durham_py",
                "recojets_durham_pz",
                "recojets_durham_e",
                "recojetconstituents_durham",
                "particlejets_durham_px",
                "particlejets_durham_py",
                "particlejets_durham_pz",
                "particlejets_durham_e",
                "particlejetconstituents_durham",
                "partonjets_durham_px",
                "partonjets_durham_py",
                "partonjets_durham_pz",
                "partonjets_durham_e",
                "partonjetconstituents_durham",

                "recojets_durhamE0_px",
                "recojets_durhamE0_py",
                "recojets_durhamE0_pz",
                "recojets_durhamE0_e",
                "recojetconstituents_durhamE0",
                "particlejets_durhamE0_px",
                "particlejets_durhamE0_py",
                "particlejets_durhamE0_pz",
                "particlejets_durhamE0_e",
                "particlejetconstituents_durhamE0",
                "recojets_durhamp_px",
                "recojets_durhamp_py",
                "recojets_durhamp_pz",
                "recojets_durhamp_e",
                "recojetconstituents_durhamp",
                "particlejets_durhamp_px",
                "particlejets_durhamp_py",
                "particlejets_durhamp_pz",
                "particlejets_durhamp_e",
                "particlejetconstituents_durhamp",

                "recojets_ee_antikt_px",
                "recojets_ee_antikt_py",
                "recojets_ee_antikt_pz",
                "recojets_ee_antikt_e",
                "recojetconstituents_ee_antikt",
                "particlejets_ee_antikt_px",
                "particlejets_ee_antikt_py",
                "particlejets_ee_antikt_pz",
                "particlejets_ee_antikt_e",
                "particlejetconstituents_ee_antikt",
                "partonjets_ee_antikt_px",
                "partonjets_ee_antikt_py",
                "partonjets_ee_antikt_pz",
                "partonjets_ee_antikt_e",
                "partonjetconstituents_ee_antikt",

                "recojets_cambridge_px",
                "recojets_cambridge_py",
                "recojets_cambridge_pz",
                "recojets_cambridge_e",
                "recojetconstituents_cambridge",
                "particlejets_cambridge_px",
                "particlejets_cambridge_py",
                "particlejets_cambridge_pz",
                "particlejets_cambridge_e",
                "particlejetconstituents_cambridge",
                "partonjets_cambridge_px",
                "partonjets_cambridge_py",
                "partonjets_cambridge_pz",
                "partonjets_cambridge_e",
                "partonjetconstituents_cambridge",

                "recojets_jade_px",
                "recojets_jade_py",
                "recojets_jade_pz",
                "recojets_jade_e",
                "recojetconstituents_jade",
                "particlejets_jade_px",
                "particlejets_jade_py",
                "particlejets_jade_pz",
                "particlejets_jade_e",
                "particlejetconstituents_jade",
                "partonjets_jade_px",
                "partonjets_jade_py",
                "partonjets_jade_pz",
                "partonjets_jade_e",
                "partonjetconstituents_jade",

                "recojets_valencia_px",
                "recojets_valencia_py",
                "recojets_valencia_pz",
                "recojets_valencia_e",
                "recojetconstituents_valencia",
                "particlejets_valencia_px",
                "particlejets_valencia_py",
                "particlejets_valencia_pz",
                "particlejets_valencia_e",
                "particlejetconstituents_valencia",
                "partonjets_valencia_px",
                "partonjets_valencia_py",
                "partonjets_valencia_pz",
                "partonjets_valencia_e",
                "partonjetconstituents_valencia",

                "Significance_reco",
                "Significance_particle",
                "Significance_parton",
                "recojets_durham_btag",
                "recojets_durhamE0_btag",
                "recojets_durhamp_btag",
                "recojets_durham_dmin",
                ]:
            branchList.push_back(branchName)
        df2.Snapshot("events", self.outname, branchList)

# example call for standalone file
# python FCCeeAnalyses/Z_Zbb_Flavor/dataframe/analysis.py /eos/experiment/fcc/ee/generation/DelphesEvents/fcc_tmp/p8_ee_Ztautau_ecm91/events_012154460.root

if __name__ == "__main__":

    if len(sys.argv)==1:
        print ("usage:")
        print ("python ",sys.argv[0]," file.root")
        sys.exit(3)
    infile = sys.argv[1]
    #outDir = '/eos/user/j/jutornda/FCCee/'+sys.argv[0].split('/')[1]+'/'
    outDir = 'outputs/semileptonic/'
    import os
    os.system("mkdir -p {}".format(outDir))
    outfile = outDir+infile.split('/')[-1]
    ncpus = 1
    analysis = analysis(infile, outfile, ncpus)
    analysis.run()

    tf = ROOT.TFile(infile)
    entries = tf.events.GetEntries()
    p = ROOT.TParameter(int)( "eventsProcessed", entries)
    outf=ROOT.TFile(outfile,"UPDATE")
    p.Write()
