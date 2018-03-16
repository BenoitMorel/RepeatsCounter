#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "pll.h"

using namespace std;

pll_msa_t *loadMSA(const string &msaFile) {
  pll_phylip_t * fd = pll_phylip_open(msaFile.c_str(), pll_map_phylip);
  if (!fd) {
    cerr << "[ERROR] Cannot open " << msaFile << endl;
    return 0;
  }
  pll_msa_t *msa = pll_phylip_parse_sequential(fd);
  if (!msa) {
    cerr << "[ERROR] Cannot parse " << msaFile << endl;
    return 0;
  }
  pll_phylip_close(fd);
  return msa;
}

void partitionMSA(const pll_msa_t *msa, 
    const string &partFile,
    vector<pll_msa_t *> &partitionnedMSAs)
{
  ifstream is(partFile);
  if (!is) {
    cerr << "[ERROR] Cannot parse " << partFile << endl;
    return;
  }
  while (is) {
    string temp;
    is >> temp;
    if (!is) 
      break;
    is >> temp; is >> temp;
    int first, last;
    is >> first;
    is.get();
    is >> last;
    first -= 1;
    pll_msa_t *submsa = (pll_msa_t *)malloc(sizeof(pll_msa_t));
    submsa->count = msa->count;
    submsa->length = last - first;
    submsa->sequence = (char **)malloc(sizeof(char*) * submsa->count);
    submsa->label = msa->label; // todobenoit do this correctly to avoid double free
    for (int i = 0; i < submsa->count; ++i) {
      submsa->sequence[i] = (char*)malloc(sizeof(char) * submsa->length);
      memcpy(submsa->sequence[i], msa->sequence[first], submsa->length);
    }
    partitionnedMSAs.push_back(submsa);
  }
}

void printHelp()
{
  cout << "./converter msa partition newick outputfile" << endl;
}

int main(int argc, char ** argv)
{
  if (argc != 5) {
    cout << "Syntax error" << endl;
    printHelp();
    exit(1);
  }
  int i = 1;
  string msaFile = argv[i++];
  string partFile = argv[i++];
  string newickFile = argv[i++];
  string outputFile = argv[i++];

  pll_msa_t *fullMSA = loadMSA(msaFile);
  vector<pll_msa_t *> partitionnedMSAs;
  partitionMSA(fullMSA, partFile, partitionnedMSAs);  

  cout << "Successfuly wrote output into " << outputFile << endl;
  return 0;
}
