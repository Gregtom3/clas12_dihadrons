import ROOT
import sys

filename = sys.argv[1]

# Open the TFile and get the TTree "EventTree"
file = ROOT.TFile(filename)
tree = file.Get("EventTree")

# Get the list of TBranches
branchList = tree.GetListOfBranches()

# Loop over the TBranches and print out the ones with the desired prefixes
for branch in branchList:
    name = branch.GetName()
    if name.startswith("gbt_") or name.startswith("rf_") or name.startswith("xgb_"):
        print(name)