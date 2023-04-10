#include "ParseBinYAML.C"
#include "injectDihadronAsym.C"

int test(){
    const char * binFile="/work/clas12/users/gmat/clas12/clas12_dihadrons/utils/binning_files/Binning.yaml";
    auto binStructs = get_structs(binFile);
    auto binStruct = binStructs[0];
    std::vector<TF2> ALUs_sig = get_ALUs_sig(binStruct);
    cout << ALUs_sig.at(0).GetTitle() << endl;
    cout << ALUs_sig.at(1).Eval(0.148099,0) << endl;
    cout << ALUs_sig.at(2).Eval(0.148099,0) << endl;
    return 0;
}