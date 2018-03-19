#include "RepeatsCounter.hpp"

#include <iostream>
#include <vector>
#include <fstream>
#include <map>

using namespace std;


using Partition = vector<vector<int> >;
using PartitionsMap = map<string, Partition>;

void parseRepeatsFile(const string &repeatsFile,
    PartitionsMap &partitionsMap)
{
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
    Partition partition(nodesNumber);
    for (int n = 0; n < nodesNumber; ++n) {
      partition[n].resize(nodesNumber);
      for (int s = 0; s < sitesNumber; ++s) {
        is >> partition[n][s];
      }
    }
    partitionsMap[partitionName] = partition;
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


  return 0;
}

