struct Block {
    std::string title;
    std::vector<std::string> var;
    std::vector<float> vmin;
    std::vector<float> vmax;
};

void writeBlockToFile(const std::string& filepath, const Block block) {
    std::ofstream file(filepath);

    if (!file) {
        std::cerr << "Error: could not open file " << filepath << " for writing" << std::endl;
        return;
    }

    // Write the block title
    file << block.title << ":" << std::endl;

    // Write each variable line
    for (size_t i = 0; i < block.var.size(); i++) {
        file << "  " << block.var[i] << ": " << block.vmin[i] << " " << block.vmax[i] << std::endl;
    }

    // Add an empty line after each block
    file << std::endl;

}

std::vector<Block> readYaml(const std::string& filepath) {
    std::ifstream file(filepath);
    std::vector<Block> blocks;

    if (!file) {
        std::cerr << "Error: could not open file " << filepath << std::endl;
        return blocks;
    }

    std::string line;
    Block currentBlock;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Ignore empty or commented lines
        }

        if (line[0] == ' ' && line != "") {
            // This is a variable line
            std::istringstream iss(line);
            std::string varName;
            float vmin, vmax;

            iss >> varName >> vmin >> vmax;

            currentBlock.var.push_back(varName.substr(0,varName.length()-1));
            currentBlock.vmin.push_back(vmin);
            currentBlock.vmax.push_back(vmax);

        } else {
            // This is a title line
            if (!currentBlock.title.empty()) {
                // We've finished reading a block, add it to the vector
                blocks.push_back(currentBlock);
            }

            currentBlock.title = line.substr(0, line.length() - 1);
            currentBlock.var.clear();
            currentBlock.vmin.clear();
            currentBlock.vmax.clear();
        }
    }

    // Add the last block we read
    if (!currentBlock.title.empty()) {
        blocks.push_back(currentBlock);
    }

    return blocks;
}