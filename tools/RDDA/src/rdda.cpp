#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>

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
  double getWeight() const { return _weight; }
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

class CPUAssignment {
public:
  CPUAssignment() {}
  void assign(const string &partitionName,
      int start,
      int size) {
    _partitions[partitionName] = pair<int, int>(start, size);
  }

private:
  map<string, pair<int, int> > _partitions;

};

bool comparePartitions(const Partition &a, 
    const Partition &b) { return (a.getWeight()  < b.getWeight()); }

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
    cout << partition.getWeight() << endl;
  }
}

void cyclicLoadBalance(vector<Partition> &partitions,
    int cores,
    vector<CPUAssignment> assignments)
{
  cout << " to implement " << endl;
}

void writeAssignment(vector<CPUAssignment> &assignments, 
    const string &outputFilename)
{
  cout << " to implement " << endl;
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
  vector<CPUAssignment> assignments;
  parseRepeatsFile(repeatsFilename, partitions);
  computeWeights(partitions);
  sort(partitions.begin(), partitions.end(), comparePartitions);
  cyclicLoadBalance(partitions, cores, assignments);

  return 0;
}



