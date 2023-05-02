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
        boundaries.push_back(parseBoundary(expr));
        str_boundaries.push_back(expr);
    }

    void insertChild(const std::shared_ptr<BinRegion>& child) {
        // Child adopts parent boundary functions
        for (const auto& boundary : boundaries) {
            child->boundaries.push_back(boundary);
        }
        // Child adopts parent boundary strings
        for (const auto& boundary : str_boundaries) {
            child->str_boundaries.push_back(boundary);
        }
        
        children.push_back(child);
    }

    bool isInside(const std::map<string, double>& coords) const {
        for (const auto& boundary : boundaries) {
            if (!boundary(coords)) {
                return false;
            }
        }
        return true;
    }

    void print() const {
        std::cout << "(  ";
        for (const auto& boundary : str_boundaries) {
            std::cout << boundary << "  ";
        }
        std::cout << ")" << std::endl;
    }
    
    void print(std::ostream& out) const {
        for (const auto& boundary : str_boundaries) {
            out << boundary << "\n";
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
    using BoundaryFunction = std::function<bool(const std::map<std::string, double>&)>;
    std::vector<BoundaryFunction> boundaries;
    std::vector<string> str_boundaries;
    std::vector<std::shared_ptr<BinRegion>> children;

    static BoundaryFunction parseBoundary(const std::string& expr) {
        auto replaceVariables = [](const std::string& expr, const std::map<std::string, double>& coords) {
            std::string replacedExpr = expr;
            for (const auto& coord : coords) {
                std::string var = coord.first;
                std::string value = std::to_string(coord.second);
                size_t pos = 0;
                while ((pos = replacedExpr.find(var, pos)) != std::string::npos) {
                    replacedExpr.replace(pos, var.length(), value);
                    pos += value.length();
                }
            }
            return replacedExpr;
        };

        return [replaceVariables, expr](const std::map<std::string, double>& coords) {
            std::string replacedExpr = replaceVariables(expr, coords);
            TFormula formula("my_formula", replacedExpr.c_str());
            // Evaluate the formula with the given coordinates
            return formula.Eval(0) > 0.0;
        };        
    }

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
};
