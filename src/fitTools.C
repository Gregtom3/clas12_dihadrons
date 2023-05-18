// Precomputed Legendre polynomials
unordered_map<string, string> precomputed_polynomials = {
    {"0_0", "1"},
    {"1_0", "cos(@th[])"},
    {"1_1", "sin(@th[])"},
    {"2_0", "1/2*(3*pow(cos(@th[]),2)-1)"},
    {"2_1", "sin(2*@th[])"}, // ignore factor 3/2
    {"2_2", "pow(sin(@th[]),2)"} // ignore factor 3
};

string legendre_polynomial(int l, int m) {
    if (m > l || l > 2) {
        // Handle invalid input
        return "Invalid input";
    }
    if(m<0){
        // Ignore prefactor in front of m<0 Associated Legendre Polynomials
        m=abs(m);
    }
    // Lookup the precomputed polynomial
    string key = to_string(l) + "_" + to_string(m);
    if (precomputed_polynomials.count(key) == 0) {
        // Handle missing precomputed polynomial
        return "Polynomial not found";
    } else {
        string polynomial = precomputed_polynomials[key];
        return polynomial;
    }
}




float get_polarization(std::string version){
  float _polarization=1.0;
  if(version=="Fall2018_RGA_inbending"){
    _polarization=0.8592;
  }
  else if(version=="Fall2018_RGA_outbending"){
    _polarization=0.8922;
  }
  else if(version=="Spring2019_RGA_inbending"){
    _polarization=0.8453;
  }
  else{
      cout << "Unknown version " << version << "...Setting polarization to 1..." << endl;
      _polarization=1.0;
  } 
  return _polarization;
}

pair<vector<string>, vector<string>> get_azi_modulations(int L, std::string version, std::string hel="hel"){
    
  float _polarization=get_polarization(version);
    
  vector<string> char_vec;
  vector<string> str_vec;
  for (int l = 0; l <= L; l++)
    {
      for (int m = 1; m <= l; m++)
        {
            
	  string str = string(Form("%f",_polarization))+"*@"+hel+"[]*sin(" + to_string(m) + "*@phi_h[]-" + to_string(m) +"*@phi_R0[])";
	  str_vec.push_back(str);
        }
      for (int m = -l; m <= l; m++)
        {
	  string str = string(Form("%f",_polarization))+"*@"+hel+"[]*sin(" + to_string(1-m) +"*@phi_h[]+" + to_string(m) +"*@phi_R0[])";
	  str_vec.push_back(str);
        }
    }
    
  // Remove duplicate entries
  sort(str_vec.begin(), str_vec.end());
  str_vec.erase(unique(str_vec.begin(), str_vec.end()), str_vec.end());
    
  // Add the mod line to the front
  for(unsigned int i = 0 ; i < str_vec.size(); i++){
    str_vec.at(i)="mod" + to_string(i) + "=" + str_vec.at(i);
  }
    
  int cidx=0;
  for (char c = 'A'; cidx<str_vec.size(); c++) 
    {
      string str = "";
      str += c;
      char_vec.push_back(str);
      cidx++;
    }
    
  return make_pair(char_vec, str_vec); 
    
}


pair<vector<string>, vector<string>> get_PW_modulations(int L, std::string version, std::string hel="hel"){
    
  float _polarization=get_polarization(version);
    
  vector<string> char_vec;
  vector<string> str_vec;
  for (int l = 0; l <= L; l++)
    {
      for (int m = 1; m <= l; m++)
        {
            
	  string str = legendre_polynomial(l, m) + "*" + string(Form("%f",_polarization))+"*@"+hel+"[]*sin(" + to_string(m) + "*@phi_h[]-" + to_string(m) +"*@phi_R0[])";
	  str_vec.push_back(str);
        }
      for (int m = -l; m <= l; m++)
        {
	  string str = legendre_polynomial(l, m) + "*" + string(Form("%f",_polarization))+"*@"+hel+"[]*sin(" + to_string(1-m) +"*@phi_h[]+" + to_string(m) +"*@phi_R0[])";
	  str_vec.push_back(str);
        }
    }
    
  // Remove duplicate entries
  sort(str_vec.begin(), str_vec.end());
  str_vec.erase(unique(str_vec.begin(), str_vec.end()), str_vec.end());
    
  // Add the mod line to the front
  for(unsigned int i = 0 ; i < str_vec.size(); i++){
    str_vec.at(i)="mod" + to_string(i) + "=" + str_vec.at(i);
  }
    
  int cidx=0;
  for (char c = 'A'; cidx<str_vec.size(); c++) 
    {
      string str = "";
      str += c;
      char_vec.push_back(str);
      cidx++;
    }
    
  return make_pair(char_vec, str_vec); 
    
}







pair<vector<string>, vector<string>> get_2h_modulations(int L, std::string version, std::string hel="hel"){

  float _polarization=get_polarization(version);
    
  vector<string> char_vec;
  vector<string> str_vec;
  for (int l = 1; l <= L; l++)
    {
      string str = string(Form("%f",_polarization))+"*@"+hel+"[]*sin(" + to_string(l) + "*@delta_phi_h[])";
      str_vec.push_back(str);
    }
    
  // Remove duplicate entries
  sort(str_vec.begin(), str_vec.end());
  str_vec.erase(unique(str_vec.begin(), str_vec.end()), str_vec.end());
    
  // Add the mod line to the front
  for(unsigned int i = 0 ; i < str_vec.size(); i++){
    str_vec.at(i)="mod" + to_string(i) + "=" + str_vec.at(i);
  }
    
  int cidx=0;
  for (char c = 'A'; cidx<str_vec.size(); c++) 
    {
      string str = "";
      str += c;
      char_vec.push_back(str);
      cidx++;
    }
    
  return make_pair(char_vec, str_vec); 
    
}






void process_azi_FM(FitManager &FM, std::string version, std::string hel="hel"){  

  auto mods = get_azi_modulations(2,version,hel);
  vector<string> char_vec = mods.first;
  vector<string> str_vec = mods.second;
    
  ///////////////////////////////Load Variables
  FM.SetUp().LoadVariable("phi_h[-3.14159265,3.14159265]");
  FM.SetUp().LoadVariable("phi_R0[-3.14159265,3.14159265]");
  FM.SetUp().LoadCategory(Form("%s[Polp=1,Polm=-1]",hel.c_str()));
  FM.SetUp().SetIDBranchName("fggID");
  ////////////////////////////////Cuts

  ///////////////////////////////Load parameters
  for (string cc: char_vec)
    FM.SetUp().LoadParameter(Form("%s[0.0,-1,1]",cc.c_str()));
  ///////////////////////////////Load formulas
  for (string ss: str_vec)
    FM.SetUp().LoadFormula(ss);
  ///////////////////////////////Load PDF
  string pdf = Form("RooComponentsPDF::AziFit(1,{phi_h,phi_R0,%s},=",hel.c_str());
  for(unsigned int i = 0 ; i < str_vec.size() ; i++){
    pdf+=char_vec.at(i);
    pdf+=";mod";
    pdf+=to_string(i);
    if(i==str_vec.size()-1)
      pdf+=")";
    else
      pdf+=":";
  }
  FM.SetUp().FactoryPDF(pdf);
  FM.SetUp().LoadSpeciesPDF("AziFit",1);
    
  return ;
}

void process_PW_FM(FitManager &FM, std::string version, std::string hel="hel"){  

  auto mods = get_PW_modulations(2,version,hel);
  vector<string> char_vec = mods.first;
  vector<string> str_vec = mods.second;
    
  ///////////////////////////////Load Variables
  FM.SetUp().LoadVariable("phi_h[-3.14159265,3.14159265]");
  FM.SetUp().LoadVariable("phi_R0[-3.14159265,3.14159265]");
  FM.SetUp().LoadVariable("th[-3.14159265,3.14159265]");
  FM.SetUp().LoadCategory(Form("%s[Polp=1,Polm=-1]",hel.c_str()));
  FM.SetUp().SetIDBranchName("fggID");
  ////////////////////////////////Cuts

  ///////////////////////////////Load parameters
  for (string cc: char_vec)
    FM.SetUp().LoadParameter(Form("%s[0.0,-1,1]",cc.c_str()));
  ///////////////////////////////Load formulas
  for (string ss: str_vec)
    FM.SetUp().LoadFormula(ss);
  ///////////////////////////////Load PDF
  string pdf = Form("RooComponentsPDF::PWFit(1,{phi_h,phi_R0,th,%s},=",hel.c_str());
  for(unsigned int i = 0 ; i < str_vec.size() ; i++){
    pdf+=char_vec.at(i);
    pdf+=";mod";
    pdf+=to_string(i);
    if(i==str_vec.size()-1)
      pdf+=")";
    else
      pdf+=":";
  }
  FM.SetUp().FactoryPDF(pdf);
  FM.SetUp().LoadSpeciesPDF("PWFit",1);
    
  return ;
}

void process_2h_FM(FitManager &FM, std::string version, std::string hel="hel"){  
  auto mods = get_2h_modulations(2,version,hel);
  vector<string> char_vec = mods.first;
  vector<string> str_vec = mods.second;
    
  ///////////////////////////////Load Variables
  FM.SetUp().LoadVariable("delta_phi_h[-3.14159265,3.14159265]");
  FM.SetUp().LoadCategory(Form("%s[Polp=1,Polm=-1]",hel.c_str()));
  FM.SetUp().SetIDBranchName("fggID");
  ////////////////////////////////Cuts

  ///////////////////////////////Load parameters
  for (string cc: char_vec)
    FM.SetUp().LoadParameter(Form("%s[0.0,-1,1]",cc.c_str()));
  ///////////////////////////////Load formulas
  for (string ss: str_vec)
    FM.SetUp().LoadFormula(ss);
  ///////////////////////////////Load PDF
  string pdf = Form("RooComponentsPDF::h2fit(1,{delta_phi_h,%s},=",hel.c_str());
  for(unsigned int i = 0 ; i < str_vec.size() ; i++){
    pdf+=char_vec.at(i);
    pdf+=";mod";
    pdf+=to_string(i);
    if(i==str_vec.size()-1)
      pdf+=")";
    else
      pdf+=":";
  }
  FM.SetUp().FactoryPDF(pdf);
  FM.SetUp().LoadSpeciesPDF("h2fit",1);
  return ;
}
