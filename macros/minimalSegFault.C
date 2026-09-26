int minimalSegFault(){
    //const char * hipoFile = "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/nSidis_005032.hipo";
    const char * hipoFile = "/cache/clas12/rg-c/production/summer22/pass1/10.5gev/NH3/dst/train/sidisdvcs/sidisdvcs_016164.hipo";
    clas12root::HipoChain _chain;
    clas12::clas12reader *_config_c12{nullptr};

    _chain.Add(hipoFile);
    _config_c12=_chain.GetC12Reader();
    _config_c12->applyQA("pass1");
    _config_c12->db()->qadb_addQARequirement("OkForAsymmetry");
    //auto &_c12=_chain.C12ref();
    //_c12->db()->qadb_requireOkForAsymmetry(true);  
    return 0;
}