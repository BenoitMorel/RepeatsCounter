#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include <stack>
#include <cmath>

using namespace std;

class Partition {
public:
  Partition() {}
  Partition(const string &name, int nodesNumber, int sitesNumber) :
    _name(name),
    _repeats(nodesNumber, vector<int>(sitesNumber)),
    _weight(0) {
      if (!nodesNumber || !sitesNumber) {
        cerr << "Error: empty partition" << endl;
      } 
    }

  int getIdentifier (int node, int site) const {return _repeats[node][site];}
  int &getIdentifier (int node, int site) {return _repeats[node][site];}
  int getNodesNumber() const {return _repeats.size();}
  int getSitesNumber() const {return _repeats[0].size();}
  const string &getName() const {return _name;}
  double getTotalWeight() const { return _weight; }
  double getPerSiteWeight() const { return _weight / getSitesNumber(); }
  void computeWeight() {
    _weight = 0.0;
    for (auto &clv: _repeats) {
      int maxId = 0;
      for (auto id: clv) {
        maxId = max(id, maxId);
      }
      _weight += (maxId + 1);
    }  
  }
private:
  string _name;
  vector<vector<int> > _repeats;
  double _weight;
};

class CoreAssignment {
public:
  CoreAssignment() : _index(0), _isFull(false) {}
  CoreAssignment(int index) : _index(index), _isFull(false) {}
  int assign(const string &partitionName,
      int start,
      int size,
      double partitionPerSiteWeight) {
    _partitions[partitionName] = pair<int, int>(start, size);
    _weight += double(size) * partitionPerSiteWeight;
  }

  double getWeight() const { return _weight; }
  bool isFull() const { return _isFull; } 
  void setFull() { _isFull = true; }
  friend ostream& operator<<(ostream& os, const CoreAssignment& assignment);  

private:
  map<string, pair<int, int> > _partitions;
  double _weight;
  int _index;
  bool _isFull;
};
  
ostream& operator<<(ostream& os, const CoreAssignment& assignment)  
{
  os << "Core" << assignment._index << " " << assignment._partitions.size() << endl; 
  for (auto p: assignment._partitions) {
    os << p.first << " "  << p.second.second;
    for (unsigned int i = 0; i < p.second.second; ++i) {
      os << " " << i + p.second.first;
    }
    os << endl;
  }
}

bool comparePartitions(const Partition &a, 
    const Partition &b) { return (a.getTotalWeight()  < b.getTotalWeight()); }

void parseRepeatsFile(const string &repeatsFile,
    vector<Partition>  &partitions)
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
    partitions.push_back(Partition(partitionName, nodesNumber, sitesNumber));
    Partition &partition = partitions[partitions.size() - 1];
    for (int n = 0; n < nodesNumber; ++n) {
      for (int s = 0; s < sitesNumber; ++s) {
        is >> partition.getIdentifier(n, s);
      }
    }
  }
  cout << "ok!" << endl;
}

void computeWeights(vector<Partition> &partitions)
{
  for (auto &partition: partitions) {
    partition.computeWeight();
  }
}

void cyclicLoadBalance(vector<Partition> &partitions,
    int cores,
    vector<CoreAssignment> &assignments)
{
  assignments.clear();
  for (unsigned int i = 0; i < cores; ++i) {
    assignments.push_back(CoreAssignment(i));
  }
  // compute the limit weight
  double limit = 0;
  for (auto &partition: partitions) {
    limit += partition.getTotalWeight() / cores;
  }
  cout << "Weight limit: " << limit << endl;
  // cyclic assignment without split
  int currentPartition = 0;
  int currentCore = 0;
  for (;currentPartition < partitions.size(); ++currentPartition) {
    Partition &partition = partitions[currentPartition];
    CoreAssignment &core = assignments[currentCore];
    if (partition.getTotalWeight() + core.getWeight() > limit) {
      break;
    }
    core.assign(partition.getName(), 0, partition.getSitesNumber(), partition.getPerSiteWeight());
    currentCore = (currentCore + 1) % cores;
  }

  // cyclic split
  for (; currentPartition < partitions.size(); ++currentPartition) {
    Partition &partition = partitions[currentPartition];
    int partitionOffset = 0;
    while (partitionOffset < partition.getSitesNumber()) {
      CoreAssignment &core = assignments[currentCore];
      int remainingSites = partition.getSitesNumber() - partitionOffset;
      double remainingWeight = remainingSites * partition.getPerSiteWeight();
      double freeWeight = limit - core.getWeight();
      int sitesToAssign = 0;
      cout << "remaining weight " << remainingWeight << endl;
      cout << "free weight " << freeWeight << endl;
      if (remainingWeight < freeWeight) {
        // assign the whole partition
        sitesToAssign = remainingSites;
      } else {
        // fill the whole core
        sitesToAssign = ceil(freeWeight / partition.getPerSiteWeight());
        core.setFull();
      }
      core.assign(partition.getName(), partitionOffset, sitesToAssign, partition.getPerSiteWeight());
      partitionOffset += sitesToAssign;
      if (assignments[currentCore].isFull()) { // would be better to keep a container of not full cores
        currentCore = (currentCore + 1) % cores;
      }
    }
  }
}

void rddaLoadBalance(vector<Partition> &partitions,
    int cores,
    vector<CoreAssignment> &assignments)
{
  assignments.clear();
  for (unsigned int i = 0; i < cores; ++i) {
    assignments.push_back(CoreAssignment(i));
  }
  // compute the limit weight
  double limit = 0;
  for (auto &partition: partitions) {
    limit += partition.getTotalWeight() / cores;
  }
  cout << "Weight limit: " << limit << endl;
  // cyclic assignment without split
  int currentPartition = 0;
  int currentCore = 0;
  for (;currentPartition < partitions.size(); ++currentPartition) {
    Partition &partition = partitions[currentPartition];
    CoreAssignment &core = assignments[currentCore];
    if (partition.getTotalWeight() + core.getWeight() > limit) {
      break;
    }
    core.assign(partition.getName(), 0, partition.getSitesNumber(), partition.getPerSiteWeight());
    currentCore = (currentCore + 1) % cores;
  }
  stack<int> qlowCores;
  stack<int> qhighCores;
  stack<int> &qlow = qlowCores;
  stack<int> &qhigh = qhighCores;

  for (unsigned int c = 0; c < currentCore; ++c) {
    qhigh.push(c);
  }
  for (unsigned int c = currentCore; c < cores; ++c) {
    qlow.push(c);
  }

 
  for (; currentPartition < partitions.size(); ++currentPartition) {
    Partition &partition = partitions[currentPartition];
    int partitionRemainingSites = partition.getSitesNumber();
    while (partitionRemainingSites > 0) {
      int partitionOffset = partition.getSitesNumber() - partitionRemainingSites;
      CoreAssignment &core = assignments[qhigh.top()];
      int sitesToAssign = ceil( (limit - core.getWeight()) / partition.getPerSiteWeight());
      if (partitionRemainingSites > sitesToAssign) {
        core.assign(partition.getName(), partitionOffset, sitesToAssign, partition.getPerSiteWeight()); 
        partitionRemainingSites -= sitesToAssign;
        qhigh.pop();
      } else {
        core = assignments[qlow.top()];
        sitesToAssign = ceil( (limit - core.getWeight()) / partition.getPerSiteWeight());
        if (partitionRemainingSites > sitesToAssign) {
          core.assign(partition.getName(), partitionOffset, sitesToAssign, partition.getPerSiteWeight()); 
          partitionRemainingSites -= sitesToAssign;
          qlow.pop();
        } else {
          core = assignments[qlow.top()];
          sitesToAssign = partitionRemainingSites;
          core.assign(partition.getName(), partitionOffset, sitesToAssign, partition.getPerSiteWeight()); 
          partitionRemainingSites -= sitesToAssign;
          qhigh.push(qlow.top()); 
          qlow.pop();
        }
      }
      if (!qlow.size()) { // trick: when qlow is empty, always fill qhigh instead
        qlow = qhigh;
      }
    
    
    }
  }

    /*
  // cyclic split
  for (; currentPartition < partitions.size(); ++currentPartition) {
    Partition &partition = partitions[currentPartition];
    int partitionOffset = 0;
    while (partitionOffset < partition.getSitesNumber()) {
      CoreAssignment &core = assignments[currentCore];
      int remainingSites = partition.getSitesNumber() - partitionOffset;
      double remainingWeight = remainingSites * partition.getPerSiteWeight();
      double freeWeight = limit - core.getWeight();
      int sitesToAssign = 0;
      cout << "remaining weight " << remainingWeight << endl;
      cout << "free weight " << freeWeight << endl;
      if (remainingWeight < freeWeight) {
        // assign the whole partition
        sitesToAssign = remainingSites;
      } else {
        // fill the whole core
        sitesToAssign = ceil(freeWeight / partition.getPerSiteWeight());
        core.setFull();
      }
      core.assign(partition.getName(), partitionOffset, sitesToAssign, partition.getPerSiteWeight());
      partitionOffset += sitesToAssign;
      if (assignments[currentCore].isFull()) { // would be better to keep a container of not full cores
        currentCore = (currentCore + 1) % cores;
      }
    }
  }
  */
}

void writeAssignment(vector<CoreAssignment> &assignments, 
    const string &outputFilename)
{
  ofstream os(outputFilename);
  os <<  assignments.size() << endl;
  int index = 0;
  for (auto assignment: assignments) {
    os << assignment;
  }
}

void printHelp()
{
  cout << "Syntax: rdda repeats_file cores_number output_file" << endl;
}

int main(int argc, char ** argv)
{
  if (argc != 4) {
    printHelp();
    return 1;
  }
  int i = 1;
  string repeatsFilename = argv[i++];
  int cores = atoi(argv[i++]);
  string outputFilename = argv[i++];
 
  vector<Partition> partitions;
  vector<CoreAssignment> assignments;
  parseRepeatsFile(repeatsFilename, partitions);
  computeWeights(partitions);
  sort(partitions.begin(), partitions.end(), comparePartitions);
  cyclicLoadBalance(partitions, cores, assignments);
  writeAssignment(assignments, outputFilename);
  return 0;
}



