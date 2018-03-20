#include "RepeatsCounter.hpp"

#include <iostream>
#include <vector>
#include <fstream>
#include <map>

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

void printHelp() {
  cout << "Syntax:" << endl;
  cout << "./RepeatsCounter repeats_file distribution_file" << endl;
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
  vector<Core> cores;
  parseDistributionFile(distributionFile, cores);
  bool consistent = checkConsistency(partitionsMap, cores);
  if (!consistent) {
    cerr << "Inconsistent data distribution, please check that each site is assigned once and only once" << endl;
    return 1;
  }
  return 0;
}

