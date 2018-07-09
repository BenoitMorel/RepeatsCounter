#include "RepeatsCounter.hpp"

#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <set>
using namespace std;

class Partition {
public:
  Partition() {}
  Partition(const string &name, int nodesNumber, int sitesNumber) :
    _name(name),
    _repeats(nodesNumber, vector<int>(sitesNumber)) {
      if (!nodesNumber || !sitesNumber) {
        cerr << "Error: empty partition" << endl;
      } 
    }

  int getIdentifier (int node, int site) const {return _repeats[node][site];}
  int &getIdentifier (int node, int site) {return _repeats[node][site];}
  int getNodesNumber() const {return _repeats.size();}
  int getSitesNumber() const {return _repeats[0].size();}
  const string &getName() const {return _name;}
private:
  string _name;
  vector<vector<int> > _repeats;
};

using PartitionsMap = map<string, Partition>;

struct SubPartition {
  string partitionName;
  vector<int> sites;
};

struct Core {
  Core() {}
  string name;
  vector<SubPartition> subPartitions;
};

bool checkConsistency(const PartitionsMap &partitionsMap, 
    const vector<Core> &cores) 
{
  map<string, vector<char> > consistencyVectorsMap;
  // prepare consistency vectors
  for (auto &partitionCell: partitionsMap) {
    auto &partition = partitionCell.second;
    consistencyVectorsMap[partition.getName()] = vector<char>(partition.getSitesNumber(), 1);
  }

  // fill consistency vectors
  for (auto &core: cores) {
    for (auto &subPartition: core.subPartitions) {
      auto &consistencyVector = consistencyVectorsMap[subPartition.partitionName];
      for (int site: subPartition.sites) {
        consistencyVector[site]--;
      }
    }
  }

  // check consistency vectors
  for (auto &partitionCell: partitionsMap) {
    auto &partition = partitionCell.second;
    if (consistencyVectorsMap[partition.getName()] != vector<char>(partition.getSitesNumber(), 0)) {
      cerr << "Unconsistency with partition " << partition.getName() << endl;
      return false;
    }
  }
  return true;
}

int getCoreCost(const PartitionsMap &partitionsMap, const Core &core)
{
  int cost = 0;
  for (auto &subPartition: core.subPartitions) {
    const Partition &partition = partitionsMap.at(subPartition.partitionName);
    for (int n = 0; n < partition.getNodesNumber(); n++) {
      set<int> identifiers;
      for (auto site: subPartition.sites) {
        identifiers.insert(partition.getIdentifier(n, site));
      }
      cost += identifiers.size();
    }
  }
  return cost;
}

int getCoreSites(Core core) {
  int sites = 0;
  for (auto &subpartition: core.subPartitions) {
    sites += subpartition.sites.size();
  }
  return sites;
}

void analyseDistribution(const PartitionsMap &partitionsMap,
    const vector<Core> &cores,
    vector<int> &costs,
    vector<int> &partitionNumbers,
    int &worstCost,
    int &totalCost
    ) 
{
  costs.clear();
  partitionNumbers.clear();
  worstCost= 0;
  totalCost = 0;
  for (auto &core: cores) {
    int cost = getCoreCost(partitionsMap, core);
    worstCost = std::max(worstCost, cost);
    totalCost += cost;
    costs.push_back(cost);
    partitionNumbers.push_back(core.subPartitions.size());
    /*
      cout << core.name << " RCC: "  << cost 
         << "\t sites: " << getCoreSites(core)
         << "\t partitions: " << core.subPartitions.size()  << endl;
         */
  }
  /*
  cout << "Worst RCC: \t\t" << worst << endl;
  cout << "Worst RCC * cores: \t" << worst * cores.size() << endl;
  cout << "Sum of RCCs: \t \t" << totalCost << endl;
  */
}

void parseRepeatsFile(const string &repeatsFile,
    PartitionsMap &partitionsMap)
{
  cout << "Parsing repeats file " << repeatsFile << "..." << endl;
  ifstream is(repeatsFile);
  if (!is) {
    cerr << "[Error] Cannot find repeats file " << repeatsFile << endl;
    return;
  }
  int partitionsNumber = 0;
  int nodesNumber = 0; 

  is >> partitionsNumber;
  is >> nodesNumber;

  for (int p = 0; p < partitionsNumber; ++p) {
    string partitionName;
    int sitesNumber = 0;
    is >> partitionName;
    is >> sitesNumber;
    partitionsMap[partitionName] = Partition(partitionName, nodesNumber, sitesNumber);
    Partition &partition = partitionsMap[partitionName];
    for (int n = 0; n < nodesNumber; ++n) {
      for (int s = 0; s < sitesNumber; ++s) {
        is >> partition.getIdentifier(n, s);
      }
    }
  }
}

void parseDistributionFile(const string &distributionFile, 
    vector<Core> &cores) 
{
  cout << "Parsing distribution file " << distributionFile << "..." << endl;
  ifstream is(distributionFile);
  if (!is) {
    cerr << "[Error} Cannot find distribution file " << distributionFile << endl;
    return;
  }
  int coresNumber = 0;
  is >> coresNumber;
  cores.resize(coresNumber);
  for (auto &core: cores) {
    is >> core.name;
    int partitionsNumber = 0;
    is >> partitionsNumber;
    core.subPartitions.resize(partitionsNumber);
    for (auto &subPartition: core.subPartitions) {
      is >> subPartition.partitionName;
      int sites = 0;
      is >> sites;
      subPartition.sites.resize(sites);
      for (auto &site: subPartition.sites) {
        is >> site;
      }
    }
  }
}

void generateOneCoreDistribution(PartitionsMap &partitions, 
    vector<Core> &cores)
{
  cores.push_back(Core());
  Core &core = cores[0];
  core.name = "UniqueCore";
  for (auto &partition: partitions) {
    core.subPartitions.push_back(SubPartition());
    SubPartition &sub = core.subPartitions[core.subPartitions.size() - 1];
    sub.partitionName = partition.second.getName();
    for (int i = 0; i < partition.second.getSitesNumber(); ++i) {
      sub.sites.push_back(i);
    }
  }
}

void printHelp() {
  cout << "Syntax:" << endl;
  cout << "./RepeatsCounter repeats_file distribution_file" << endl;
}

void analyse(PartitionsMap &partitionsMap, vector<Core> &cores,
    vector<int> &costs,
    vector<int> &partitionNumbers,
    int &worstCore,
    int &totalCost)
{
  bool consistent = checkConsistency(partitionsMap, cores);
  if (!consistent) {
    cerr << "Inconsistent data distribution, please check that each site is assigned once and only once" << endl;
    return;
  }
  analyseDistribution(partitionsMap, cores, costs, partitionNumbers, worstCore, totalCost);
}

void printTitle(const string &msg) 
{
  cout << "***********************************" << endl;
  cout << "** \t" << msg << "\t **" << endl;
  cout << "***********************************" << endl;
  
}

int main(int argc, char **argv)
{
  if (argc != 3) {
    cerr << "[Error] Invalid number of arguments" << endl;
    printHelp();
    return 1;
  }

  int i = 1;
  string repeatsFile = argv[i++];
  string distributionFile = argv[i++];
  
  PartitionsMap partitionsMap;
  parseRepeatsFile(repeatsFile, partitionsMap);
  vector<Core> uniqueCore;
  vector<Core> cores;
  generateOneCoreDistribution(partitionsMap, uniqueCore);
  parseDistributionFile(distributionFile, cores);
  //printTitle("Analyse with one single core");
  //
  vector<int> singleCosts;
  vector<int> singlePartitions;
  vector<int> costs;
  vector<int> partitions;
  int singleWorstCore;
  int singleTotalCost;
  int worstCore;
  int totalCost;
  printTitle("Analysis results");
  analyse(partitionsMap, uniqueCore, singleCosts, singlePartitions, singleWorstCore, singleTotalCost);
  analyse(partitionsMap, cores, costs, partitions, worstCore, totalCost);

  int maxPartitions = 0;
  
  for (unsigned int i = 0; i < cores.size(); ++i) {
    auto &core = cores[i];
    maxPartitions = max(maxPartitions, (int)core.subPartitions.size());
    cout << core.name << " RCC: "  << costs[i] 
         << "\t sites: " << getCoreSites(core) 
         << "\t partitions: " << core.subPartitions.size()  << endl;
  }
  cout << endl;
  cout << "Worst RCC lower bound:\t" << singleTotalCost / cores.size() << " (total RCC before split / number of cores)" << endl;
  cout << "Worst RCC: \t\t" << worstCore << endl;
  cout << "Worst RCC * cores: \t" << worstCore * cores.size() << endl;
  cout << "Sum of RCCs: \t \t" << totalCost << endl;
  cout << "Total repeats loss:\t" << totalCost - singleTotalCost << endl;
  cout << "Worst partitions:\t" << maxPartitions << endl;
  return 0;
}

