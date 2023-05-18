#include <iostream>
#include <vector>
#include <functional>
#include <string>
#include <stdexcept>
#include <map>
#include <memory>
#include <atomic>

class BinRegion {
public:
    BinRegion() : id(nextId++) {}

    void addBoundary(const std::string& expr) {
        TFormula formula(("formula_" + std::to_string(id) + "_" + std::to_string(boundaryFormulas.size())).c_str(), expr.c_str());
        boundaryFormulas.push_back(formula);
    }

    void insertChild(const std::shared_ptr<BinRegion>& child) {
        // Child adopts parent boundaries
        for (auto& formula : boundaryFormulas) {
            child->boundaryFormulas.push_back(formula);
        }
        children.push_back(child);
    }

    bool isInside(const std::map<std::string, double>& coords) const {
        for (TFormula formula : boundaryFormulas) {
            for (const auto& coord : coords) {
                formula.SetParameter(coord.first.c_str(), coord.second);
            }
            if (formula.Eval(0) <= 0.0) {
                return false;
            }
        }
        return true;
    }

    void print() const {
        std::cout << "(  ";
        for (const auto& boundary : boundaryFormulas) {
            std::cout << boundary.GetExpFormula() << "  ";
        }
        std::cout << ")" << std::endl;
    }
    
    void print(std::ostream& out) const {
        for (const auto& boundary : boundaryFormulas) {
            out << boundary.GetExpFormula() << "\n";
        }
    }
    
    int getId() const {
        return id;
    }

    std::shared_ptr<BinRegion> getChild(size_t index) const {
        if (index < children.size()) {
            return children[index];
        }
        return nullptr;
    }
    
    bool hasChildren() const {
        return !children.empty();
    }

    size_t numChildren() const {
        return children.size();
    }

    static std::atomic<int> nextId;

private:
    int id;
    std::vector<TFormula> boundaryFormulas;
    std::vector<std::shared_ptr<BinRegion>> children;
};

std::atomic<int> BinRegion::nextId(0);

std::vector<std::shared_ptr<BinRegion>> createEvenlySpacedBins(string dimension, int numBins, double minValue, double maxValue) {
    std::vector<std::shared_ptr<BinRegion>> bins;

    double binSize = (maxValue - minValue) / numBins;
    for (int i = 0; i < numBins; ++i) {
        auto bin = std::make_shared<BinRegion>();

        std::string lowerBoundary = dimension + (i == 0 ? ">" : ">=") + std::to_string(minValue + i * binSize);
        std::string upperBoundary = dimension + (i == numBins - 1 ? "<" : "<=") + std::to_string(minValue + (i + 1) * binSize);

        bin->addBoundary(lowerBoundary);
        bin->addBoundary(upperBoundary);
        bins.push_back(bin);
    }

    return bins;
}


std::vector<std::shared_ptr<BinRegion>> createCustomBins(string dimension, const std::vector<float>& binEdges) {
    std::vector<std::shared_ptr<BinRegion>> bins;

    for (size_t i = 0; i < binEdges.size() - 1; ++i) {
        auto bin = std::make_shared<BinRegion>();

        std::string lowerBoundary = dimension + (i == 0 ? ">" : ">=") + std::to_string(binEdges[i]);
        std::string upperBoundary = dimension + (i == binEdges.size() - 2 ? "<" : "<=") + std::to_string(binEdges[i + 1]);

        bin->addBoundary(lowerBoundary);
        bin->addBoundary(upperBoundary);
        bins.push_back(bin);
    }

    return bins;
}

class BinMap {
public:
    void setRoot(const std::shared_ptr<BinRegion>& rootRegion) {
        root = rootRegion;
    }

    std::shared_ptr<BinRegion> getRoot() const {
        return root;
    }

    void serialize(std::ostream& out) const {
        serializeRegion(root, out);
    }
    
    int findLowestLevelChild(const std::map<std::string, double>& coords) const {
        return findLowestLevelChild(coords, root);
    }
    
private:
    std::shared_ptr<BinRegion> root;

    void serializeRegion(const std::shared_ptr<BinRegion>& region, std::ostream& out) const {
        if (region->hasChildren()) {
            for (size_t i = 0; i < region->numChildren(); ++i) {
                serializeRegion(region->getChild(i), out);
            }
        } else {
            out << "unique id: " << region->getId() << std::endl;
            region->print(out);
        }
    }
    
    int findLowestLevelChild(const std::map<std::string, double>& coords, const std::shared_ptr<BinRegion>& region) const {
        if (region->hasChildren()) {
            int foundChildId = -1;
            for (size_t i = 0; i < region->numChildren(); ++i) {
                int childId = findLowestLevelChild(coords, region->getChild(i));
                if (childId != -1) {
                    if (foundChildId == -1) {
                        foundChildId = childId;
                    } else {
                        std::cerr << "Error: Data falls into more than one lowest level child boundary." << std::endl;
                        return -1;
                    }
                }
            }
            return foundChildId;
        } else {
            if (region->isInside(coords)) {
                return region->getId();
            } else {
                return -1;
            }
        }
    }
};
