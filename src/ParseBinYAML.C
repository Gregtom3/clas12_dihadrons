struct YAMLbinstruct {
    std::string name = "";
    int numDimensions = 0;
    int nBins = 1;
    std::vector<std::string> dimensionNames;
    std::vector<std::vector<double>> binEdges;
    std::string injectName = "";
    std::vector<std::string> injectsigFuncs;
    std::vector<std::string> injectbgFuncs;
};

// Function prototypes
std::vector<YAMLbinstruct> get_structs(const char* inyaml);
int findIndexOfFile(const char* file, const std::vector<YAMLbinstruct>& binStructs);
YAMLbinstruct get_struct(const char* inyaml, const char* file);
void ParseBinYaml(const char* input_data_wildcard, const char* output_datadir, const char* input_yaml, int isMC, int binNum);
void serializeYAMLbinstruct(const YAMLbinstruct& binstruct, const std::string& filename);


void printYAMLbinstruct(const YAMLbinstruct& binstruct) {
    std::cout << "name: " << binstruct.name << std::endl;
    std::cout << "numDimensions: " << binstruct.numDimensions << std::endl;
    std::cout << "nBins: " << binstruct.nBins << std::endl;
    std::cout << "dimensionNames: ";
    for (const auto& name : binstruct.dimensionNames) {
        std::cout << name << " ";
    }
    std::cout << std::endl;
    std::cout << "binEdges: ";
    for (const auto& edges : binstruct.binEdges) {
        for (const auto& edge : edges) {
            std::cout << edge << " ";
        }
        std::cout << "| ";
    }
    std::cout << std::endl;
    std::cout << "injectName: " << binstruct.injectName << std::endl;
    std::cout << "injectsigFuncs: ";
    for (const auto& func : binstruct.injectsigFuncs) {
        std::cout << func << " ";
    }
    std::cout << std::endl;
    std::cout << "injectbgFuncs: ";
    for (const auto& func : binstruct.injectbgFuncs) {
        std::cout << func << " ";
    }
    std::cout << std::endl;
}

std::vector<std::string> parseList(const std::string& line) {
    std::vector<std::string> result;
    std::string content = line.substr(line.find("[") + 1, line.find("]") - line.find("[") - 1);
    content.erase(std::remove_if(content.begin(), content.end(), ::isspace), content.end());
    while (content.find(",") != std::string::npos) {
        result.push_back(content.substr(0, content.find(",")));
        content = content.substr(content.find(",") + 1);
    }
    result.push_back(content);
    return result;
}

std::vector<double> parseVectorDouble(const std::string& line) {
    std::vector<double> result;
    std::string content = line.substr(line.find("[") + 1, line.find("]") - line.find("[") - 1);
    while (content.find(",") != std::string::npos) {
        result.push_back(std::stod(content.substr(0, content.find(","))));
        content = content.substr(content.find(",") + 1);
    }
    result.push_back(std::stod(content));
    return result;
}

std::vector<YAMLbinstruct> get_structs(const char* inyaml) {
    std::ifstream yamlFile(inyaml);
    std::vector<YAMLbinstruct> YAMLbinstructs;

    // Read line by line
    std::string line;
    YAMLbinstruct currentStructure;

    while (std::getline(yamlFile, line)) {
        if (line.find("name") != std::string::npos) {
            if (currentStructure.name != "") {
                for (const auto& vect : currentStructure.binEdges) {
                    currentStructure.nBins *= (vect.size() - 1);
                }
                YAMLbinstructs.push_back(currentStructure);
                currentStructure = {};
            }
            currentStructure.name = line.substr(line.find(":") + 2);
        } else if (line.find("numDimensions") != std::string::npos) {
            currentStructure.numDimensions = std::stoi(line.substr(line.find(":") + 2));
        } else if (line.find("dimensionNames") != std::string::npos) {
            currentStructure.dimensionNames = parseList(line);
        } else if (line.find("binEdges") != std::string::npos) {
            while (std::getline(yamlFile, line) && line.find("-") != std::string::npos) {
                currentStructure.binEdges.push_back(parseVectorDouble(line));
            }
        } 
        if (line.find("injectLabel") != std::string::npos) {
            currentStructure.injectName = line.substr(line.find(":") + 2);
        } else if (line.find("inject_sigfuncs") != std::string::npos) {
            while (std::getline(yamlFile, line) && line.find("-") != std::string::npos) {
                currentStructure.injectsigFuncs.push_back(line.substr(line.find("\"") + 1, line.find_last_of("\"") - line.find("\"") - 1));
            }
        } 
        if (line.find("inject_bgfuncs") != std::string::npos) {
            while (std::getline(yamlFile, line) && line.find("-") != std::string::npos) {
                currentStructure.injectbgFuncs.push_back(line.substr(line.find("\"") + 1, line.find_last_of("\"") - line.find("\"") - 1));
            }
        }
    }

    // Include this line to save the last structure
    if (currentStructure.name != "") {
                for (const auto& vect : currentStructure.binEdges) {
                    currentStructure.nBins *= (vect.size() - 1);
                }
                YAMLbinstructs.push_back(currentStructure);
                currentStructure = {};
    }
    yamlFile.close();
    
    return YAMLbinstructs;
}

void serializeYAMLbinstruct(const YAMLbinstruct& binstruct, const std::string& filename) {
    std::ofstream outfile(filename);
    outfile << "name: " << binstruct.name << std::endl;
    outfile << "numDimensions: " << binstruct.numDimensions << std::endl;
    outfile << "nBins: " << binstruct.nBins << std::endl;
    outfile << "dimensionNames: ";
    for (const auto& name : binstruct.dimensionNames) {
        outfile << name << " ";
    }
    outfile << std::endl;
    outfile << "binEdges: ";
    for (const auto& edges : binstruct.binEdges) {
        for (const auto& edge : edges) {
            outfile << edge << " ";
        }
        outfile << "| ";
    }
    outfile << std::endl;
    outfile << "injectName: " << binstruct.injectName << std::endl;
    outfile << "injectsigFuncs: ";
    for (const auto& func : binstruct.injectsigFuncs) {
        outfile << func << " ";
    }
    outfile << std::endl;
    outfile << "injectbgFuncs: ";
    for (const auto& func : binstruct.injectbgFuncs) {
        outfile << func << " ";
    }
    outfile << std::endl;
    outfile.close();
}


int findIndexOfFile(const char* file, const std::vector<YAMLbinstruct>& binStructs) {
    int index = -1;
    std::string fileName = file;
    int pos = fileName.find_last_of("/\\");
    if (fileName.find("nSidis") != std::string::npos) {
        fileName = fileName.substr(pos + 8);
    } else {
        fileName = fileName.substr(pos + 4);
    }
    for (int i = 0; i < binStructs.size(); i++) {
        if (binStructs[i].name == fileName) {
            index = i;
            break;
        }
    }
    return index;
}

YAMLbinstruct get_struct(const char* inyaml, const char* file) { // Get the binstruct corresponding to the specifically binned file
    auto bs = get_structs(inyaml);
    return bs.at(findIndexOfFile(file, bs));
}
