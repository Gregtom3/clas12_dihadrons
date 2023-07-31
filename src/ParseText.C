
//Function to get pid_h1 and pid_h2 from filename
void getPIDs(std::string filename, int& pid_h1, int& pid_h2, std::string &particleNames) {

    // Extract the string between the second and third underscores
    std::size_t first = filename.find("/data/");
    std::size_t second = filename.find('/', first+6);
    particleNames = filename.substr(first+6,second-first-6);

    // Transform the particle names to lower case
    std::transform(particleNames.begin(), particleNames.end(), particleNames.begin(), ::tolower);
    // Assign the appropriate PID values to pid_h1 and pid_h2
    if (particleNames == "piplus_pi0") {
        pid_h1 = 211;
        pid_h2 = 111;
    } else if (particleNames == "piminus_pi0") {
        pid_h1 = -211;
        pid_h2 = 111;
    } else if (particleNames == "pi0_pi0") {
        pid_h1 = 111;
        pid_h2 = 111;
    } else if (particleNames == "piplus_piminus"){
        pid_h1 = 211;
        pid_h2 = -211;
    } else if (particleNames == "piplus_piplus"){
        pid_h1 = 211;
        pid_h2 = 211;
    } else if (particleNames == "piminus_piminus"){
        pid_h1 = -211;
        pid_h2 = -211;
    } else {
        // Invalid input string, set PIDs to zero
        pid_h1 = 0;
        pid_h2 = 0;
    }
}

// Input --> "/path/to/file/<version>_merged.root"
// Output --> "<version>"
std::string getVersion(const char *c =""){
    std::string s(c);
    std::size_t lastSlash = s.find_last_of('/');
    std::string filename = s.substr(lastSlash + 1);
    std::string prefix = filename.substr(0, filename.find("_merged.root"));
    return prefix;
}

// Input --> "/path/to/file/<version>_merged_cuts.root"
// Output --> "<version>"
std::string getVersion2(const char *c =""){
    std::string s(c);
    std::size_t lastSlash = s.find_last_of('/');
    std::string filename = s.substr(lastSlash + 1);
    std::string prefix = filename.substr(0, filename.find("_merged_cuts.root"));
    return prefix;
}

// Input --> "/path/to/file/nSidis_RGA_5032.root"
// Output --> "5032"
int extractRunNumberFromDataFile(const char *c) {
    std::string filepath(c);
    
    // find the last _ in the string
    size_t underscorePos = filepath.rfind("_");

    if (underscorePos == std::string::npos) {
        std::cerr << "Invalid filepath: " << filepath << std::endl;
        return -1;
    }

    // find the .root extension in the string
    size_t rootPos = filepath.rfind(".root");

    if (rootPos == std::string::npos) {
        std::cerr << "Invalid filepath: " << filepath << std::endl;
        return -1;
    }

    // make sure the number is left of the .root
    if (underscorePos > rootPos) {
        std::cerr << "Invalid filepath: " << filepath << std::endl;
        return -1;
    }

    // Start of the number
    size_t startPos = underscorePos + 1;


    // get the substring that represents the number
    std::string numberStr = filepath.substr(startPos, rootPos - startPos);

    // convert the number string to an integer and return
    return std::stoi(numberStr);
}


std::string get_dir_from_binstruct_idx(YAMLbinstruct binStruct, int binnum){
    std::string dir="";
    
    int numDimensions = binStruct.numDimensions;
    std::vector<std::vector<double>> binEdges = binStruct.binEdges;

    std::vector<int> indices(numDimensions, 0);
    
    int M = 0;
    while (true) {
        // Clear the values vector
        dir="";
        bool skip = false;
        // Append the minimum and maximum bin edges for each dimension to the values vector
        for (int i = 0; i < numDimensions; i++) {
            int index = indices[i];
            double minVal = binEdges[i][index];
            double maxVal = (index+1 < binEdges[i].size()) ? binEdges[i][index+1] : binEdges[i][index];

            // Check if the minimum and maximum bin edges are equal
            if (minVal == maxVal) {
                dir="";
                i = numDimensions;
                M+=-1;
                // Skip this dimension if the minimum and maximum bin edges are equal
                continue;
            }
            
            dir+=binStruct.dimensionNames.at(i)+"_"+std::to_string(minVal)+"_"+std::to_string(maxVal);
            if(i+1!=numDimensions){
                dir+="_";
            }
        }
        
        if(M==binnum){
            return dir;
        }
        else{
            M+=1;
        }

        // Increment the indices
        int dimIndex = numDimensions - 1;
        while (dimIndex >= 0) {
            if (indices[dimIndex] < binEdges[dimIndex].size() - 1) {
                indices[dimIndex]++;
                break;
            }
            indices[dimIndex] = 0;
            dimIndex--;
        }

        // Check if we have printed out all bin combinations
        if (dimIndex < 0) {
            break;
        }
    }
    return "";
}



std::string get_cut_from_binstruct_idx(YAMLbinstruct binStruct, int binnum){
    std::string cut="";
    
    int numDimensions = binStruct.numDimensions;
    std::vector<std::vector<double>> binEdges = binStruct.binEdges;

    std::vector<int> indices(numDimensions, 0);
    
    int M = 0;
    while (true) {
        // Clear the values vector
        cut="";
        bool skip = false;
        // Append the minimum and maximum bin edges for each dimension to the values vector
        for (int i = 0; i < numDimensions; i++) {
            int index = indices[i];
            double minVal = binEdges[i][index];
            double maxVal = (index+1 < binEdges[i].size()) ? binEdges[i][index+1] : binEdges[i][index];

            // Check if the minimum and maximum bin edges are equal
            if (minVal == maxVal) {
                cut="";
                i = numDimensions;
                M+=-1;
                // Skip this dimension if the minimum and maximum bin edges are equal
                continue;
            }
            
            cut+=binStruct.dimensionNames.at(i)+">"+std::to_string(minVal)+"&&"+binStruct.dimensionNames.at(i)+"<"+std::to_string(maxVal);
            if(i+1!=numDimensions){
                cut+="&&";
            }
        }
        
        if(M==binnum){
            return cut;
        }
        else{
            M+=1;
        }

        // Increment the indices
        int dimIndex = numDimensions - 1;
        while (dimIndex >= 0) {
            if (indices[dimIndex] < binEdges[dimIndex].size() - 1) {
                indices[dimIndex]++;
                break;
            }
            indices[dimIndex] = 0;
            dimIndex--;
        }

        // Check if we have printed out all bin combinations
        if (dimIndex < 0) {
            break;
        }
    }
    return "";
}