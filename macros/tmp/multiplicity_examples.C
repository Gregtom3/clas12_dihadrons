#include "../src/MultiplicityBinning.C"

int multiplicity_examples() {
    auto rootRegion = std::make_shared<BinRegion>();
    rootRegion->addBoundary("x+y>3");
    rootRegion->addBoundary("x>0");
    rootRegion->addBoundary("x<1");
    rootRegion->addBoundary("y>-3");

    auto childRegion1 = std::make_shared<BinRegion>();
    childRegion1->addBoundary("z>0");
    childRegion1->addBoundary("z<0.5");

    auto childRegion2 = std::make_shared<BinRegion>();
    childRegion2->addBoundary("z>=0.5");
    childRegion2->addBoundary("z<=1.0");

    rootRegion->insertChild(childRegion1);
    rootRegion->insertChild(childRegion2);
    
    BinMap binMap;
    binMap.setRoot(rootRegion);

    // Test if a point is inside a region
    std::map<std::string, double> point = {{"x", 0.5}, {"y", 3}, {"z", 0.3}};
    std::cout << "Is the point (0.5, 0.5, 0.3) inside rootRegion? "
              << (binMap.getRoot()->isInside(point) ? "Yes" : "No") << std::endl;
    std::cout << "Is the point (0.5, 0.5, 0.3) inside childRegion1? "
              << (childRegion1->isInside(point) ? "Yes" : "No") << std::endl;
    std::cout << "Is the point (0.5, 0.5, 0.3) inside childRegion2? "
              << (childRegion2->isInside(point) ? "Yes" : "No") << std::endl;

//     auto rootRegion = std::make_shared<BinRegion>();
    
//     std::vector<float> wBinEdges = {0.0, 0.2, 0.5, 1.0};
//     std::vector<std::shared_ptr<BinRegion>> wBins = createCustomBins("w", wBinEdges);
//     for (const auto& wBin : wBins) {
//         rootRegion->insertChild(wBin);
//     }
    
//     rootRegion->getChild(0)->print();
    
//     BinMap binMap;
//     binMap.setRoot(rootRegion);
    
    std::ofstream outFile("bin_map_output.txt");
    if (outFile.is_open()) {
        binMap.serialize(outFile);
        outFile.close();
    } else {
        std::cerr << "Unable to open output file." << std::endl;
    }
    return 0;
}