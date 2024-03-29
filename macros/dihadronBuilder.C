#include "../src/Constants.h"
#include "../src/HipoBankInterface.C"
#include "../src/CLAS12Analysis.C"
#include "../src/Kinematics.C"
#include "../src/ParseBinYAML.C"
#include "../src/ParseText.C"

// This program reads in a root file containing data from CLAS12 collisions and the corresponding particles produced. 
// It then uses the information to construct dihadrons (two hadrons) from the particles. 
// It reads in the necessary information from the root file (generated by hipo2tree.C) and uses it to determine the kinematics of the dihadron and fill a new TTree with the information. 
// The TTree contains the kinematics of each dihadron and is saved back into the file

// The program uses the input_file name to determine what hadrons to build. If pi0's are built, then the "weight_branch" tells the program which machine learning model
// is to be used to save photon classification values. The "weight_branch" in the "EventTree" is created by the program "machine_learning/photonID/predict.py"

int dihadronBuilder(const char *input_file="hipo2tree.root",
                    const char *weight_branch="none"){
    

    // Declare pid_h1 and pid_h2
    int pid_h1=0;
    int pid_h2=0;
    std::string particleNames ="";
    // Determine the pids from the file (see function)
    getPIDs(string(input_file),pid_h1,pid_h2,particleNames);
    // Read the TFile
    TFile *f = new TFile(input_file,"UPDATE");
    // Read the TTree
    TTree *EventTree = (TTree*)f->Get("EventTree");
    // Declare important branches
    //declare all necessary variables
    double x, Q2, W, Pol,y;
    double truex, trueQ2, trueW, truey;
    int hel,run,A,_evnum;
    int Nmax=100;
    double px[Nmax], py[Nmax], pz[Nmax], E[Nmax], vz[Nmax], chi2[Nmax], theta[Nmax], eta[Nmax], phi[Nmax];
    double truepx[Nmax] , truepy[Nmax] , truepz[Nmax], trueE[Nmax], truetheta[Nmax], trueeta[Nmax], truephi[Nmax];
    double weight[Nmax];
    int parentID[Nmax],parentPID[Nmax],parentparentID[Nmax],parentparentPID[Nmax];
    int is_CFR[Nmax];
    int pid[Nmax], truepid[Nmax];
    //link the TBranches to the variables
    EventTree->SetBranchAddress("A",&A);
    EventTree->SetBranchAddress("evnum",&_evnum);
    EventTree->SetBranchAddress("run",&run);
    EventTree->SetBranchAddress("Pol",&Pol);
    EventTree->SetBranchAddress("hel",&hel);
    EventTree->SetBranchAddress("trueparentid", &parentID);
    EventTree->SetBranchAddress("trueparentpid", &parentPID);
    EventTree->SetBranchAddress("trueparentparentid", &parentparentID);
    EventTree->SetBranchAddress("trueparentparentpid", &parentparentPID);
    EventTree->SetBranchAddress("x",&x);
    EventTree->SetBranchAddress("Q2",&Q2);
    EventTree->SetBranchAddress("W",&W); 
    EventTree->SetBranchAddress("y",&y); 
    EventTree->SetBranchAddress("Nmax",&Nmax);
    EventTree->SetBranchAddress("px",px);
    EventTree->SetBranchAddress("py",py);
    EventTree->SetBranchAddress("pz",pz);
    EventTree->SetBranchAddress("E",E);
    EventTree->SetBranchAddress("vz",vz);
    EventTree->SetBranchAddress("chi2",chi2);
    EventTree->SetBranchAddress("pid",pid);
    EventTree->SetBranchAddress("theta",theta);
    EventTree->SetBranchAddress("eta",eta);
    EventTree->SetBranchAddress("phi",phi);
    EventTree->SetBranchAddress("truex",&truex);
    EventTree->SetBranchAddress("truey",&truey);
    EventTree->SetBranchAddress("trueQ2",&trueQ2);
    EventTree->SetBranchAddress("trueW",&trueW);
    EventTree->SetBranchAddress("trueE",&trueE);
    EventTree->SetBranchAddress("truepx",&truepx);
    EventTree->SetBranchAddress("truepy",&truepy);
    EventTree->SetBranchAddress("truepz",&truepz);
    EventTree->SetBranchAddress("truetheta",&truetheta);
    EventTree->SetBranchAddress("trueeta",&trueeta);
    EventTree->SetBranchAddress("truephi",&truephi);
    EventTree->SetBranchAddress("truepid",truepid);
    EventTree->SetBranchAddress("is_CFR",is_CFR);
    if((pid_h1==111||pid_h2==111)&&std::string(weight_branch)!="none")
        EventTree->SetBranchAddress(weight_branch,weight);
    
    //Create the new TTree
    TString treename="dihadron";
    //If dihadron tree already exists, remove it
    if (f->Get(treename)) f->Delete("dihadron*;*");
    TTree *outtree = new TTree(treename,"Dihadron-by-Dihadron info");
    
    double M1,M2,Mh,phi_h,phi_R0,phi_R1,th,z1,z2,xF1,xF2,z,xF,Mx,phi_h1,phi_h2,delta_phi_h,pT_1,pT_2,pT_tot,P_1,P_2,P_tot;
    double  trueM1,trueM2,trueMh,truephi_h,truephi_R0,truephi_R1,trueth,truez1,truez2,truexF1,truexF2,truez,truexF,trueMx,truephi_h1,truephi_h2,truedelta_phi_h,truepT_1,truepT_2,truepT_tot,trueP_1,trueP_2,trueP_tot;
    int truepid_e;
    double E_e, th_e, phi_e;
    int truepid_1,truepid_2,trueparentpid_1,trueparentpid_2,trueparentid_1,trueparentid_2,trueparentparentpid_1,trueparentparentpid_2,trueparentparentid_1,trueparentparentid_2, trueparentpid_11, trueparentpid_12, trueparentpid_21, trueparentpid_22;
    double E_11, E_12, E_21, E_22;
    double th_11, th_12, th_21, th_22;
    double phi_11, phi_12, phi_21, phi_22;

    int is_CFR_1, is_CFR_2;
    int MCmatch; // MCmatch --> 1 if all particles have Monte Carlo pairing
    int isGoodEventWithoutML;
    int isGoodEventWithML;
    int truepid_11, truepid_12, truepid_21, truepid_22; // For photon pairs
    double trueM12, M12; // addition of M1 M2
    double fgID=0; //uniqueID for each dihadron
    // Machine Learning output
    double p_11=-1;
    double p_12=-1;
    double p_21=-1;
    double p_22=-1;
    // Create branches
    outtree->Branch("A", &A, "A/I");
    outtree->Branch("evnum", &_evnum, "evnum/I");
    outtree->Branch("fgID", &fgID, "fgID/D");
    outtree->Branch("run", &run, "run/I");
    outtree->Branch("Pol", &Pol, "Pol/D");
    outtree->Branch("hel", &hel, "hel/I");
    outtree->Branch("MCmatch", &MCmatch, "MCmatch/I");
    outtree->Branch("isGoodEventWithoutML", &isGoodEventWithoutML, "isGoodEventWithoutML/I");
    outtree->Branch("isGoodEventWithML", &isGoodEventWithML, "isGoodEventWithML/I");
    outtree->Branch("is_CFR_1",&is_CFR_1, "is_CFR_1/I");
    outtree->Branch("is_CFR_2",&is_CFR_2, "is_CFR_2/I");
    outtree->Branch("x", &x, "x/D");
    outtree->Branch("Q2", &Q2, "Q2/D");
    outtree->Branch("W", &W, "W/D");
    outtree->Branch("y", &y, "y/D");
    outtree->Branch("M1", &M1, "M1/D");
    outtree->Branch("M2", &M2, "M2/D");
    outtree->Branch("M12",&M12,"M12/D");
    outtree->Branch("Mh", &Mh, "Mh/D");
    outtree->Branch("phi_h", &phi_h, "phi_h/D");
    outtree->Branch("phi_R0", &phi_R0, "phi_R0/D");
    outtree->Branch("phi_R1", &phi_R1, "phi_R1/D");
    outtree->Branch("th", &th, "th/D");
    outtree->Branch("z1", &z1, "z1/D");
    outtree->Branch("z2", &z2, "z2/D");
    outtree->Branch("xF1", &xF1, "xF1/D");
    outtree->Branch("xF2", &xF2, "xF2/D");
    outtree->Branch("z", &z, "z/D");
    outtree->Branch("xF", &xF, "xF/D");
    outtree->Branch("Mx", &Mx, "Mx/D");
    outtree->Branch("E_e", &E_e, "E_e/D");
    outtree->Branch("th_e", &th_e, "th_e/D");
    outtree->Branch("phi_e", &phi_e, "phi_e/D");
    outtree->Branch("phi_h1", &phi_h1, "phi_h1/D");
    outtree->Branch("phi_h2", &phi_h2, "phi_h2/D");
    outtree->Branch("delta_phi_h", &delta_phi_h, "delta_phi_h/D");
    outtree->Branch("pT1", &pT_1, "pT1/D");
    outtree->Branch("pT2", &pT_2, "pT2/D");
    outtree->Branch("pTtot", &pT_tot, "pTtot/D");
    outtree->Branch("P1", &P_1, "P1/D");
    outtree->Branch("P2", &P_2, "P2/D");
    outtree->Branch("Ptot", &P_tot, "Ptot/D");
    outtree->Branch("truex", &truex, "truex/D");
    outtree->Branch("trueQ2", &trueQ2, "trueQ2/D");
    outtree->Branch("trueW", &trueW, "trueW/D");
    outtree->Branch("truey", &truey, "truey/D");
    outtree->Branch("trueM1", &trueM1, "trueM1/D");
    outtree->Branch("trueM2", &trueM2, "trueM2/D");
    outtree->Branch("trueM12",&trueM12,"trueM12/D");
    outtree->Branch("trueMh", &trueMh, "trueMh/D");
    outtree->Branch("truephi_h", &truephi_h, "truephi_h/D");
    outtree->Branch("truephi_R0", &truephi_R0, "truephi_R0/D");
    outtree->Branch("truephi_R1", &truephi_R1, "truephi_R1/D");
    outtree->Branch("trueth", &trueth, "trueth/D");
    outtree->Branch("truez1", &truez1, "truez1/D");
    outtree->Branch("truez2", &truez2, "truez2/D");
    outtree->Branch("truexF1", &truexF1, "truexF1/D");
    outtree->Branch("truexF2", &truexF2, "truexF2/D");
    outtree->Branch("truez", &truez, "truez/D");
    outtree->Branch("truexF", &truexF, "truexF/D");
    outtree->Branch("trueMx", &trueMx, "trueMx/D");
    outtree->Branch("truephi_h1", &truephi_h1, "truephi_h1/D");
    outtree->Branch("truephi_h2", &truephi_h2, "truephi_h2/D");
    outtree->Branch("truedelta_phi_h", &truedelta_phi_h, "truedelta_phi_h/D");
    outtree->Branch("truepT1", &truepT_1, "truepT1/D");
    outtree->Branch("truepT2", &truepT_2, "truepT2/D");
    outtree->Branch("truepTtot", &truepT_tot, "truepTtot/D");
    outtree->Branch("trueP1", &trueP_1, "trueP1/D");
    outtree->Branch("trueP2", &trueP_2, "trueP2/D");
    outtree->Branch("truePtot", &trueP_tot, "truePtot/D");
    outtree->Branch("truepid_e",&truepid_e, "truepid_e/I");
    outtree->Branch("truepid_1", &truepid_1, "truepid_1/I");
    outtree->Branch("truepid_2", &truepid_2, "truepid_2/I");
    outtree->Branch("truepid_11", &truepid_11, "truepid_11/I");
    outtree->Branch("truepid_12", &truepid_12, "truepid_12/I");
    outtree->Branch("truepid_21", &truepid_21, "truepid_21/I");
    outtree->Branch("truepid_22", &truepid_22, "truepid_22/I");
    outtree->Branch("trueparentpid_1", &trueparentpid_1, "trueparentpid_1/I");
    outtree->Branch("trueparentpid_2", &trueparentpid_2, "trueparentpid_2/I");
    outtree->Branch("trueparentpid_11", &trueparentpid_11, "trueparentpid_11/I");
    outtree->Branch("trueparentpid_12", &trueparentpid_12, "trueparentpid_12/I");
    outtree->Branch("trueparentpid_21", &trueparentpid_21, "trueparentpid_21/I");
    outtree->Branch("trueparentpid_22", &trueparentpid_22, "trueparentpid_22/I");
    outtree->Branch("trueparentid_1", &trueparentid_1, "trueparentid_1/I");
    outtree->Branch("trueparentid_2", &trueparentid_2, "trueparentid_2/I");
    outtree->Branch("trueparentparentpid_1", &trueparentparentpid_1, "trueparentparentpid_1/I");
    outtree->Branch("trueparentparentpid_2", &trueparentparentpid_2, "trueparentparentpid_2/I");
    outtree->Branch("trueparentparentid_1", &trueparentparentid_1, "trueparentparentid_1/I");
    outtree->Branch("trueparentparentid_2", &trueparentparentid_2, "trueparentparentid_2/I");
    outtree->Branch("E_11", &E_11, "E_11/D");
    outtree->Branch("E_12", &E_12, "E_12/D");
    outtree->Branch("E_21", &E_21, "E_21/D");
    outtree->Branch("E_22", &E_22, "E_22/D");
    outtree->Branch("th_11", &th_11, "th_11/D");
    outtree->Branch("th_12", &th_12, "th_12/D");
    outtree->Branch("th_21", &th_21, "th_21/D");
    outtree->Branch("th_22", &th_22, "th_22/D");
    outtree->Branch("phi_11", &phi_11, "phi_11/D");
    outtree->Branch("phi_12", &phi_12, "phi_12/D");
    outtree->Branch("phi_21", &phi_21, "phi_21/D");
    outtree->Branch("phi_22", &phi_22, "phi_22/D");
    outtree->Branch("p_11", &p_11,"p_11/D");
    outtree->Branch("p_12", &p_12,"p_12/D");
    outtree->Branch("p_21", &p_21,"p_21/D");
    outtree->Branch("p_22", &p_22,"p_22/D");
    
    // Clone the outtree to only fill it if the cuts pass
    TTree *outtree_clone = outtree->CloneTree(-1, "fast");
    outtree_clone->SetName("dihadron_cuts");
    
    // Kinematics/CLAS12Analysis Object
    Kinematics kin;
    CLAS12Analysis clas12ana = CLAS12Analysis();
    
    // Initial particles
    TLorentzVector init_electron(0,0,0,0); // To be set one run is found
    TLorentzVector init_target(0,0,0,Mp);
    
    // for loop over all events
    int N = EventTree->GetEntries();
    std::vector<std::vector<int>> dihadron_idxs;    
    for (int ev=0; ev<N; ev++){

        if((ev+1)%100==0 || ev==N-1){
            if(ev!=N-1){
                cout << "Progress: " << ev+1 << "/" << N << "\r";
                cout.flush();
            }
            else
                cout << "Progress: " << ev+1 << "/" << N << endl ;
        }
        

        EventTree->GetEntry(ev);

        if(ev==0){
            init_electron.SetE(runBeamEnergy(run));
            init_electron.SetPz(sqrt(init_electron.E()*init_electron.E()-Me*Me));
        }
        
        //Loop over all particles in the event to find electron
        TLorentzVector electron;
        TLorentzVector trueelectron;
        TLorentzVector q; // virtual photon
        TLorentzVector trueq;
        int idx_e=-1;
        double max_e=-1;
        for (int i=0; i<Nmax; i++){
            if(pid[i]==11){
                if(E[i]>max_e){
                    idx_e=i;
                    max_e=E[i];
                }
            }
        }
        
        electron.SetPxPyPzE(px[idx_e],py[idx_e],pz[idx_e],E[idx_e]);
        trueelectron.SetPxPyPzE(truepx[idx_e],truepy[idx_e],truepz[idx_e],trueE[idx_e]);
        truepid_e=truepid[idx_e];
        E_e=electron.E();
        th_e=electron.Theta();
        phi_e=electron.Phi();
        q=init_electron-electron;
        trueq=init_electron-trueelectron;
        
        dihadron_idxs = clas12ana.dihadron_idxs(pid_h1,pid_h2,pid,Nmax);
        
        // Now loop over all dihadrons
        for(int a = 0 ; a < dihadron_idxs.size() ; a++){
            std::vector<int> dihadron_idx = dihadron_idxs.at(a);
            int i=0;
            int ii=0;
            int j=0;
            int jj=0;
            if(pid_h1==111){
                i=dihadron_idx.at(0);
                ii=dihadron_idx.at(1);
            }else{
                i=dihadron_idx.at(0);
            }
            if(pid_h2==111&&pid_h1!=111){
                j=dihadron_idx.at(1);
                jj=dihadron_idx.at(2);
            }else if(pid_h2==111&&pid_h1==111){
                j=dihadron_idx.at(2);
                jj=dihadron_idx.at(3);
            }else if(pid_h1==111){
                j=dihadron_idx.at(2);
            }else{
                j=dihadron_idx.at(1);
            }
            
            TLorentzVector h1;
            TLorentzVector trueh1;
            TLorentzVector h2;
            TLorentzVector trueh2;
            TLorentzVector dihadron;
            TLorentzVector truedihadron;
  
            if(pid_h1==111){
                p_11 = weight[i];
                p_12 = weight[ii];
		E_11 = E[i];
		E_12 = E[ii];
		th_11 = theta[i];
		th_12 = theta[ii];
		phi_11 = phi[i];
		phi_12 = phi[ii];
                h1.SetPxPyPzE(px[i]+px[ii],py[i]+py[ii],pz[i]+pz[ii],E[i]+E[ii]);
                trueh1.SetPxPyPzE(truepx[i]+truepx[ii],truepy[i]+truepy[ii],truepz[i]+truepz[ii],trueE[i]+trueE[ii]);
                truepid_11=truepid[i];
                truepid_12=truepid[i];
                trueparentpid_11=parentPID[i];
                trueparentpid_12=parentPID[ii];
                if(parentID[i]==parentID[ii] && parentID[i]!=-999){
                    trueparentid_1=parentID[i];
                    trueparentpid_1=parentPID[i];
                    trueparentparentid_1=parentparentID[i];
                    trueparentparentpid_1=parentparentPID[i];
                    is_CFR_1=is_CFR[i];
                } else {
                   trueparentid_1=-999;
                   trueparentpid_1=-999;
                   trueparentparentid_1=-999;
                   trueparentparentpid_1=-999;
                   is_CFR_1=-999;
                }
            }else{
	      E_11 = E[i];
	      E_12 = -999;
	      th_11 = theta[i];
	      th_12 = -999;
	      phi_11 = phi[i];
	      phi_12 = -999;
	      h1.SetPxPyPzE(px[i],py[i],pz[i],E[i]);
                trueh1.SetPxPyPzE(truepx[i],truepy[i],truepz[i],trueE[i]);
                truepid_1=truepid[i];
                trueparentid_1=parentID[i];
                trueparentpid_1=parentPID[i];
                trueparentparentid_1=parentparentID[i];
                trueparentparentpid_1=parentparentPID[i];
                trueparentpid_11=parentPID[i];
                trueparentpid_12=parentPID[i];
                is_CFR_1=is_CFR[i];
            }
            if(pid_h2==111){
                p_21 = weight[j];
                p_22 = weight[jj];
		E_21 = E[j];
		E_22 = E[jj];
		th_21 = theta[j];
		th_22 = theta[jj];
		phi_21 = phi[j];
		phi_22 = phi[jj];
                h2.SetPxPyPzE(px[j]+px[jj],py[j]+py[jj],pz[j]+pz[jj],E[j]+E[jj]);
                trueh2.SetPxPyPzE(truepx[j]+truepx[jj],truepy[j]+truepy[jj],truepz[j]+truepz[jj],trueE[j]+trueE[jj]);
                truepid_21=truepid[j];
                truepid_22=truepid[j];
                trueparentpid_21=parentPID[j];
                trueparentpid_22=parentPID[jj];
                if(parentID[j]==parentID[jj] && parentID[j]!=-999){
                    trueparentid_2=parentID[j];
                    trueparentpid_2=parentPID[j];
                    trueparentparentid_2=parentparentID[j];
                    trueparentparentpid_2=parentparentPID[j];
                    is_CFR_2=is_CFR[j];
                } else {
                   trueparentid_2=-999;
                   trueparentpid_2=-999;
                   trueparentparentid_2=-999;
                   trueparentparentpid_2=-999;
                   is_CFR_2=-999;
                }
            }
            else{
	      E_21 = E[j];
	      E_22 = -999;
	      th_21 = theta[j];
	      th_22 = -999;
	      phi_21 = phi[j];
	      phi_22 = -999;
	      
                h2.SetPxPyPzE(px[j],py[j],pz[j],E[j]);
                trueh2.SetPxPyPzE(truepx[j],truepy[j],truepz[j],trueE[j]);
                truepid_2=truepid[j];
                trueparentid_2=parentID[j];
                trueparentpid_2=parentPID[j];
                trueparentparentid_2=parentparentID[j];
                trueparentparentpid_2=parentparentPID[j];
                trueparentpid_21=parentPID[j];
                trueparentpid_22=parentPID[j];
                is_CFR_2=is_CFR[j];
            }
            
            if(pid_h1==pid_h2){
                z1 = kin.z(init_target,h1,q);
                z2 = kin.z(init_target,h2,q);
                if(z1<z2){
                    TLorentzVector temp = h1;
                    h2=h1;
                    h1=temp;
                    TLorentzVector truetemp = h1;
                    trueh2=trueh1;
                    trueh1=truetemp;
                    std::swap(p_11,p_21);
                    std::swap(p_12,p_22);
		    std::swap(E_11, E_21);
		    std::swap(E_12, E_22);
		    std::swap(th_11, th_21);
		    std::swap(th_12, th_22);
		    std::swap(phi_11, phi_21);
		    std::swap(phi_12, phi_22);
		    std::swap(truepid_11,truepid_21);
                    std::swap(truepid_12,truepid_22);
                    std::swap(truepid_1,truepid_2);
                    std::swap(trueparentid_1,trueparentid_2);
                    std::swap(trueparentpid_1,trueparentpid_2);
                    std::swap(trueparentparentid_1,trueparentparentid_2);
                    std::swap(trueparentparentpid_1,trueparentparentpid_2);
                    std::swap(is_CFR_1,is_CFR_2);
                }
            }
            
            // Build the dihadron
            dihadron = h1+h2;
            truedihadron = trueh1+trueh2;
            // fill results
            M1 = h1.M();
            M2 = h2.M();
            M12 = M1+M2;
            Mh = dihadron.M();
            phi_h = kin.phi_h(q,init_electron,h1,h2);
            phi_h1 = kin.phi_h(q,init_electron,h1);
            phi_h2 = kin.phi_h(q,init_electron,h2);
            delta_phi_h = phi_h1-phi_h2;
            if(delta_phi_h>PI){
                delta_phi_h-=2*PI;
            }else if(delta_phi_h<-PI){
                delta_phi_h+=2*PI;
            }
            pT_1 = kin.Pt(q,h1,init_target);
            pT_2 = kin.Pt(q,h2,init_target);
            pT_tot = kin.Pt(q,dihadron,init_target);
            P_1 = h1.P();
            P_2 = h2.P();
            P_tot = dihadron.P();
            phi_R0 = kin.phi_R(q,init_electron,h1,h2,0);
            phi_R1 = kin.phi_R(q,init_electron,h1,h2,1);
            th     = kin.com_th(h1,h2);
            xF1 = kin.xF(q,h1,init_target,W);
            xF2 = kin.xF(q,h2,init_target,W);
            xF     = kin.xF(q,dihadron,init_target,W);
            z1 = kin.z(init_target,h1,q);
            z2 = kin.z(init_target,h2,q);
            z = z1+z2;
            Mx = (init_electron+init_target-electron-dihadron).M();


            trueM1 = trueh1.M();
            trueM2 = trueh2.M();
            trueM12 = trueM1+trueM2;
            trueMh = truedihadron.M();
            truephi_h = kin.phi_h(trueq,init_electron,trueh1,trueh2);
            truephi_h1 = kin.phi_h(trueq,init_electron,trueh1);
            truephi_h2 = kin.phi_h(trueq,init_electron,trueh2);
            truedelta_phi_h = truephi_h1-truephi_h2;
            if(truedelta_phi_h>PI){
                truedelta_phi_h-=2*PI;
            }else if(truedelta_phi_h<-PI){
                truedelta_phi_h+=2*PI;
            }
            truepT_1 = kin.Pt(trueq,trueh1,init_target);
            truepT_2 = kin.Pt(trueq,trueh2,init_target);
            truepT_tot = kin.Pt(trueq,truedihadron,init_target);
            trueP_1 = trueh1.P();
            trueP_2 = trueh2.P();
            trueP_tot = truedihadron.P();
            truephi_R0 = kin.phi_R(trueq,init_electron,trueh1,trueh2,0);
            truephi_R1 = kin.phi_R(trueq,init_electron,trueh1,trueh2,1);
            trueth     = kin.com_th(trueh1,trueh2);
            truexF1 = kin.xF(trueq,trueh1,init_target,trueW);
            truexF2 = kin.xF(trueq,trueh2,init_target,trueW);
            truexF     = kin.xF(trueq,truedihadron,init_target,trueW);
            truez1 = kin.z(init_target,trueh1,trueq);
            truez2 = kin.z(init_target,trueh2,trueq);
            truez = truez1+truez2;
            trueMx = (init_electron+init_target-trueelectron-truedihadron).M();
            MCmatch=0;
            if(trueelectron.E()>0&&trueh1.E()>0&&trueh2.E()>0) MCmatch=1;
            
            // Determine if this would've been a good event without ML
            isGoodEventWithoutML=1;
            if(pid_h1==111){
                isGoodEventWithoutML*=(E[i]>0.6);
                isGoodEventWithoutML*=(E[ii]>0.6);
            }
            
            if(pid_h2==111){
                isGoodEventWithoutML*=(E[j]>0.6);
                isGoodEventWithoutML*=(E[jj]>0.6);
            }
            
            // Determine if we should fill the cloned TTree if it passes our cuts
            bool fill_clone = true;
            fill_clone*=(z<0.95);
            fill_clone*=(xF1>0&&xF2>0);
            if((pid_h1==211&&pid_h2==-211)||(pid_h1==211&&pid_h2==111)){
                fill_clone*=(Mx>1.5);
            }
            if(pid_h1==211||pid_h1==-211){
                fill_clone*=(P_1>1.25);
            }
            if(pid_h2==211||pid_h2==-211){
                fill_clone*=(P_2>1.25);
            }
            if(pid_h1==111){
                fill_clone*=(p_11>0.9&&p_12>0.9);
            }
	    // Set isGoodEventWithoutML up to this point
	    isGoodEventWithoutML*=fill_clone;
	    // Now set the fill_clone if the ML cut passes
            if(pid_h2==111){
                fill_clone*=(p_21>0.9&&p_22>0.9);
            }
	    // Set isGoodEventWithML to fill_clone
	    isGoodEventWithML = fill_clone;
	    // Fill the cloned, abrigded tree
            if(fill_clone){
                outtree_clone->Fill();
            }
	    // Fill the larger, main tree without cuts
            outtree->Fill();
            fgID++;
        } // end dihadron loop
    }// end event loop

   
    cout << "Writing Total TTree with " << outtree->GetEntries() << " entries" << endl;
    cout << "Writing Cut TTree with " << outtree_clone->GetEntries() << " entries" << endl;
    outtree->Write();
    outtree_clone->Write();
    f->Close();
    cout << "Done" << endl;
    return 0;
}


