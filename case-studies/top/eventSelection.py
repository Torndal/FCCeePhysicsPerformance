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

               .Define('EVTrest_thrust',     'Algorithms::minimize_thrust("Minuit2","Migrad")(RPrest_px, RPrest_py, RPrest_pz)')
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
               
               ################## Durham ###################
               #run jet clustering with reconstructed particles. Durham_algorithm, exclusive clustering, exactly 4 jets, sorted by E, E-scheme  
               .Define("FCCAnalysesRecoJets_durham", "JetClustering::clustering_ee_kt(2, 4, 1, 0)(reco_jets)")
               .Define("jets",            "JetClusteringUtils::get_pseudoJets(FCCAnalysesRecoJets_durham)")
               .Define("jetconstituents", "JetClusteringUtils::get_constituents(FCCAnalysesRecoJets_durham)")
               .Define("jets_dmin",       "JetClusteringUtils::get_dmerge(FCCAnalysesRecoJets_durham)")
               .Define("jets_ymin",       "JetClusteringUtils::get_ymerge(FCCAnalysesRecoJets_durham)")
               .Define("jets_px",         "JetClusteringUtils::get_px(jets)")
               .Define("jets_py",         "JetClusteringUtils::get_py(jets)")
               .Define("jets_pz",         "JetClusteringUtils::get_pz(jets)")
               .Define("jets_e",          "JetClusteringUtils::get_e(jets)")
               .Define("jets_flavour",    "JetTaggingUtils::get_flavour(jets, Particle)")
               .Define("jets_btag",       "JetTaggingUtils::get_btag(jets_flavour, 0.80)")
               .Define("jets_truebtag",       "JetTaggingUtils::get_btag(jets_flavour, 1.00)")
               
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

               .Define("jetsignificance", "SemileptonicTop::VertexSignificance(jetconstituents, RPrest, EFlowTrack_1)")

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
                
                "EVT_thrust_val",
                "EVTrest_thrust_val",
                
                #"AlgoSpher",

                "jets_px",
                "jets_py",
                "jets_pz",
                "jets_e",
                "jetconstituents",
                "jetsignificance",
                "jets_btag",
                "jets_truebtag",
                "jets_dmin",
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
    outDir = '/eos/user/j/jutornda/FCCee/'+sys.argv[1].split('/')[-2]+'/'
    #outDir = 'outputs/semileptonic/test/'
    import os
    os.system("mkdir -p {}".format(outDir))
    outfile = outDir+infile.split('/')[-1]
    ncpus = 4
    analysis = analysis(infile, outfile, ncpus)
    analysis.run()

    tf = ROOT.TFile(infile)
    entries = tf.events.GetEntries()
    p = ROOT.TParameter(int)( "eventsProcessed", entries)
    outf=ROOT.TFile(outfile,"UPDATE")
    p.Write()
