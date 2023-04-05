#ifndef Structs_h
#define Structs_h
struct FS{
  int pid_h1=0;
  int pid_h2=0;
  int num_h1=0;
  int num_h2=0;
};

FS get_FS(int pid_h1, int pid_h2){
  FS fs;
  if(pid_h1==111&&pid_h2==111){ //Pi0Pi0
    fs.pid_h1=22;
    fs.num_h1=4;
    fs.pid_h2=0;
    fs.num_h2=0;
  }
  else if(pid_h1==111&&pid_h2!=111){ // One Pi0
    fs.pid_h1=22;
    fs.num_h1=2;
    fs.pid_h2=pid_h2;
    fs.num_h2=1;
  }
  else if(pid_h1!=111&&pid_h2==111){ // One Pi0
    fs.pid_h1=pid_h1;
    fs.num_h1=1;
    fs.pid_h2=22;
    fs.num_h2=2;
  }
  else if(pid_h1==pid_h2){ // PiPi
    fs.pid_h1=pid_h1;
    fs.num_h1=2;
    fs.pid_h2=0;
    fs.num_h2=0;
  }
  else if(pid_h1==0 && pid_h2==0){
    fs.pid_h1=0;
    fs.num_h1=0;
    fs.pid_h2=0;
    fs.num_h2=0;
  }else{
    fs.pid_h1=pid_h1;
    fs.num_h1=1;
    fs.pid_h2=pid_h2;
    fs.num_h2=1;
  }
  return fs;
} 

struct part{
  // Reconstructed Info
  int pindex=-999;
  int status=0;
  double px=-999;
  double py=-999;
  double pz=-999;
  double pt=-999;
  double p=-999;
  double E=-999;
  double m=-999;
  double theta=-999;
  double eta=-999;
  double phi=-999;
  int pid=-999;
  double vx=-999;
  double vy=-999;
  double vz=-999;
  double chi2=-999;
  double beta=-999;
    
    
  // MC Lund Info
  double truepx=-999;
  double truepy=-999;
  double truepz=-999;
  double truep=-999;
  double truept=-999;
  double trueE=-999;
  double truem=-999;
  double truetheta=-999;
  double trueeta=-999;
  double truephi=-999;
  double truevx=-999;
  double truevy=-999;
  double truevz=-999;
  int is_CFR=-999;
  int truepid=-999;
  int trueparentid=-999;
  int trueparentpid=-999;
  int trueparentparentid=-999;
  int trueparentparentpid=-999;
  // Calorimeter Info
  int    pcal_sector=-999;
  double pcal_e=-999;
  double pcal_x=-999;
  double pcal_y=-999;
  double pcal_z=-999;
  double pcal_lu=-999;
  double pcal_lv=-999;
  double pcal_lw=-999;
  double pcal_m2u=-999;
  double pcal_m2v=-999;
  double pcal_m2w=-999;
    
  int    ecin_sector=-999;
  double ecin_e=-999;
  double ecin_x=-999;
  double ecin_y=-999;
  double ecin_z=-999;
  double ecin_lu=-999;
  double ecin_lv=-999;
  double ecin_lw=-999;
  double ecin_m2u=-999;
  double ecin_m2v=-999;
  double ecin_m2w=-999;

  int    ecout_sector=-999;
  double ecout_e=-999;
  double ecout_x=-999;
  double ecout_y=-999;
  double ecout_z=-999;
  double ecout_lu=-999;
  double ecout_lv=-999;
  double ecout_lw=-999;
  double ecout_m2u=-999;
  double ecout_m2v=-999;
  double ecout_m2w=-999;
    
  // Drift Chamber Info
  int sector = -999;
  double traj_x1=-999;
  double traj_y1=-999;
  double traj_z1=-999;
  double traj_x2=-999;
  double traj_y2=-999;
  double traj_z2=-999;
  double traj_x3=-999;
  double traj_y3=-999;
  double traj_z3=-999;
    
};
#endif
