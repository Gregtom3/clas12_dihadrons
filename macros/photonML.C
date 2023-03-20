
int photonML(const char * input_file = "/volatile/clas12/users/gmat/clas12analysis.sidis.data/clas12_dihadrons/projects/ana_v1/data/raw/pi0_pi0/MC_3051_0.root"){

    //Define the variables "m_g" , "m_ch" , "m_nh" 
    int m_g = 3; // Number of neighboring gammas
    int m_ch = 2; // Number of neighboring charged hadrons (protons, pions, kaons)
    int m_nh = 2; // Number of neighboring neutral hadrons (neutrons)

    //Read the TFile
    TFile *f = new TFile(input_file,"UPDATE");

    //Read the TTree
    TTree *EventTree = (TTree*)f->Get("EventTree");

    //If MLInput tree already exists, remove it
    if (f->Get("MLInput")) f->Delete("MLInput;*");

    //Create the new TTree
    TTree *MLInput = new TTree("MLInput","Nearest neighbor information");

    //Define the branches in MLInput: POI (photon-of-interest)
    //                                Nearest neighbor gammas, charged hadrons, neutral hadrons, electron
    int flag=0; // 0 --> Photon is background (no MC Match), 1 --> Photon is signal (MC Match)
    double gE=0; // Photon energy
    double gEpcal=0; // Photon pcal energy
    double gTheta=0; // Photon angle
    double gm2u=0; // Photon shower shape
    double gm2v=0; // Photon shower shape
    
    double R_gamma[m_g]; // Angular distance between calo shower centers
    double dE_gamma[m_g]; // Energy difference
    double Epcal_gamma[m_g]; // Energy deposited in the pcal
    double m2u_gamma[m_g]; // Shower shape variables
    double m2v_gamma[m_g]; // Shower shape variables
    
    double R_ch[m_ch]; // Angular distance between calo shower centers
    double dE_ch[m_ch]; // Energy difference
    double Epcal_ch[m_ch]; // Energy deposited in the pcal
    double m2u_ch[m_ch]; // Shower shape variables
    double m2v_ch[m_ch]; // Shower shape variables
    
    double R_nh[m_nh]; // Angular distance between calo shower centers
    double dE_nh[m_nh]; // Energy difference
    double Epcal_nh[m_nh]; // Energy deposited in the pcal
    double m2u_nh[m_nh]; // Shower shape variables
    double m2v_nh[m_nh]; // Shower shape variables
    
    double num_photons_0_1, num_photons_0_2, num_photons_0_35;
    
    
    // Place TTree Branches
    
    MLInput->Branch("flag",&flag,"flag/I");
    MLInput->Branch("m_g",&m_g,"m_g/I");
    MLInput->Branch("m_ch",&m_ch,"m_ch/I");
    MLInput->Branch("m_nh",&m_nh,"m_nh/I");
    
    MLInput->Branch("gE",&gE,"gE/D");
    MLInput->Branch("gEpcal",&gEpcal,"gEpcal/D");
    MLInput->Branch("gTheta",&gTheta,"gTheta/D");
    MLInput->Branch("gm2u",&gm2u,"gm2u/D");
    MLInput->Branch("gm2v",&gm2v,"gm2v/D");
    
    MLInput->Branch("R_gamma",R_gamma,"R_gamma[m_g]/D");
    MLInput->Branch("dE_gamma",dE_gamma,"dE_gamma[m_g]/D");
    MLInput->Branch("Epcal_gamma",Epcal_gamma,"Epcal_gamma[m_g]/D");
    MLInput->Branch("m2u_gamma",m2u_gamma,"m2u_gamma[m_g]/D");
    MLInput->Branch("m2v_gamma",m2v_gamma,"m2v_gamma[m_g]/D");
    
    MLInput->Branch("R_ch",R_ch,"R_ch[m_ch]/D");
    MLInput->Branch("dE_ch",dE_ch,"dE_ch[m_ch]/D");
    MLInput->Branch("Epcal_ch",Epcal_ch,"Epcal_ch[m_ch]/D");
    MLInput->Branch("m2u_ch",m2u_ch,"m2u_ch[m_ch]/D");
    MLInput->Branch("m2v_ch",m2v_ch,"m2v_ch[m_ch]/D");
    
    MLInput->Branch("R_nh",R_nh,"R_nh[m_nh]/D");
    MLInput->Branch("dE_nh",dE_nh,"dE_nh[m_nh]/D");
    MLInput->Branch("Epcal_nh",Epcal_nh,"Epcal_nh[m_nh]/D");
    MLInput->Branch("m2u_nh",m2u_nh,"m2u_nh[m_nh]/D");
    MLInput->Branch("m2v_nh",m2v_nh,"m2v_nh[m_nh]/D");
    
    MLInput->Branch("num_photons_0_1",&num_photons_0_1,"num_photons_0_1/D");
    MLInput->Branch("num_photons_0_2",&num_photons_0_2,"num_photons_0_2/D");
    MLInput->Branch("num_photons_0_35",&num_photons_0_35,"num_photons_0_35/D");
    
    //Define variables to read from EventTree
    const int kNmax = 500;
    int Nmax;
    double E[kNmax], th[kNmax], trueE[kNmax], pcal_e[kNmax], pcal_m2u[kNmax], pcal_m2v[kNmax];
    double pid[kNmax], pcal_x[kNmax],pcal_y[kNmax],pcal_z[kNmax];

    //Set the address of the branches in EventTree
    EventTree->SetBranchAddress("Nmax", &Nmax);
    EventTree->SetBranchAddress("E", E);
    EventTree->SetBranchAddress("theta", th);
    EventTree->SetBranchAddress("trueE", trueE);
    EventTree->SetBranchAddress("theta", th);
    EventTree->SetBranchAddress("pid", pid);
    EventTree->SetBranchAddress("pcal_x",pcal_x);
    EventTree->SetBranchAddress("pcal_y",pcal_y);
    EventTree->SetBranchAddress("pcal_z",pcal_z);
    EventTree->SetBranchAddress("pcal_e",pcal_e);
    EventTree->SetBranchAddress("pcal_m2u",pcal_m2u);
    EventTree->SetBranchAddress("pcal_m2v",pcal_m2v);
    
    //Loop over the events in EventTree
    for (int iEvent=0; iEvent<EventTree->GetEntries(); ++iEvent) {
      EventTree->GetEntry(iEvent);
      //Loop over the particles in the event
      for (int ipart=0; ipart<Nmax; ++ipart) {
        //Check if the particle is a photon
        if (pid[ipart] == 22) {
          //Initialize the arrays
          for (int i=0; i<m_g; ++i) {
            R_gamma[i] = 0;
            dE_gamma[i] = 0;
            Epcal_gamma[i] = 0;
            m2u_gamma[i] = 0;
            m2v_gamma[i] = 0;
          }
          for (int i=0; i<m_ch; ++i){
            R_ch[i] = 0;
            dE_ch[i] = 0;
            Epcal_ch[i] = 0;
            m2u_ch[i] = 0;
            m2v_ch[i] = 0;
          }
          for (int i=0; i<m_nh; ++i){
            R_nh[i] = 0;
            dE_nh[i] = 0;
            Epcal_nh[i] = 0;
            m2u_nh[i] = 0;
            m2v_nh[i] = 0;
          }

          //Set the number of photons within R<0.1, R<0.2, R<0.35
          num_photons_0_1 = 0;
          num_photons_0_2 = 0;
          num_photons_0_35 = 0;
          
          // Set flag
          if(trueE[ipart]>0) flag=1; // trueE == -999 if no match is found
          else flag=0;
        
          // Set vars
          gE = E[ipart];
          gEpcal = pcal_e[ipart];
          gTheta = th[ipart];
          gm2u = pcal_m2u[ipart];
          gm2v = pcal_m2v[ipart];
          //Loop over the particles in the event
          for (int jpart=0; jpart<Nmax; ++jpart) {

            //Skip the same particle
            if (ipart == jpart) continue;

            //Calculate R 
            TVector3 v_1(pcal_x[ipart],pcal_y[ipart],pcal_z[ipart]);
            TVector3 v_2(pcal_x[jpart],pcal_y[jpart],pcal_z[jpart]);
            float R = v_1.Angle(v_2);
            
            if(pid[jpart]==22){
                //Count the number of photons within R<0.1, R<0.2, R<0.35
                if (R < 0.1) num_photons_0_1++;
                if (R < 0.2) num_photons_0_2++;
                if (R < 0.35) num_photons_0_35++;
                //Find the photon's nearest photon neighbors
                //Specifically, sort by smallest R value, of which there are "m_g" in the sorted list
                //The code is a bit longer because if a smaller R is found, we must bump forward all other elements in the list
                for (int i=0; i<m_g; ++i) {
                  if (R < R_gamma[i] || R_gamma[i] == 0) {
                    int j = m_g - 1;
                    while (j > i) {
                        R_gamma[j] = R_gamma[j - 1];
                        dE_gamma[j] = dE_gamma[j - 1];
                        Epcal_gamma[j] = Epcal_gamma[j - 1];
                        m2u_gamma[j] = m2u_gamma[j - 1];
                        m2v_gamma[j] = m2v_gamma[j - 1];
                        j--;
                    }
                    R_gamma[i] = R;
                    dE_gamma[i] = E[ipart] - E[jpart];
                    Epcal_gamma[i] = pcal_e[jpart];
                    m2u_gamma[i] = pcal_m2u[jpart];
                    m2v_gamma[i] = pcal_m2v[jpart];
                    break;
                  }
                }
            }
            else if(pid[jpart]==211 || pid[jpart]==-211 || pid[jpart]==2212 || pid[jpart]==-2212 || pid[jpart]==321 || pid[jpart]==-321){
                //Find the photon's nearest charged hadrons neighbors
                for (int i=0; i<m_ch; ++i) {
                  if (R < R_ch[i] || R_ch[i] == 0) {
                    int j = m_ch - 1;
                    while (j > i) {
                        R_ch[j] = R_ch[j - 1];
                        dE_ch[j] = dE_ch[j - 1];
                        Epcal_ch[j] = Epcal_ch[j - 1];
                        m2u_ch[j] = m2u_ch[j - 1];
                        m2v_ch[j] = m2v_ch[j - 1];
                        j--;
                    }
                    R_ch[i] = R;
                    dE_ch[i] = E[ipart] - E[jpart];
                    Epcal_ch[i] = pcal_e[jpart];
                    m2u_ch[i] = pcal_m2u[jpart];
                    m2v_ch[i] = pcal_m2v[jpart];
                    break;
                  }
                }
            }
            else if(pid[jpart]==2112 || pid[jpart]==-2112){
                //Find the photon's nearest neutral hadrons neighbors
                for (int i=0; i<m_nh; ++i) {
                  if (R < R_nh[i] || R_nh[i] == 0) {
                    int j = m_nh - 1;
                    while (j > i) {
                        R_nh[j] = R_nh[j - 1];
                        dE_nh[j] = dE_nh[j - 1];
                        Epcal_nh[j] = Epcal_nh[j - 1];
                        m2u_nh[j] = m2u_nh[j - 1];
                        m2v_nh[j] = m2v_nh[j - 1];
                        j--;
                    }
                    R_nh[i] = R;
                    dE_nh[i] = E[ipart] - E[jpart];
                    Epcal_nh[i] = pcal_e[jpart];
                    m2u_nh[i] = pcal_m2u[jpart];
                    m2v_nh[i] = pcal_m2v[jpart];
                    break;
                  }
                }
            }
            else{
                continue;
            }
              
            
          }

          //Fill the MLInput TTree
          MLInput->Fill();
        }
      }
    }

    //Write the MLInput TTree to disk
    MLInput->Write();

    //Close the TFile
    f->Close();
    return 0;
}