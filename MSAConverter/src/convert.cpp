#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <search.h>
#include "pll.h"

using namespace std;
pll_utree_t *loadTree(const string &newickFile) {
  pll_utree_t * tree = pll_utree_parse_newick(newickFile.c_str());
  if (!tree) {
    cerr << "[Error] could not parse " << newickFile << endl;
  }
  hcreate(tree->tip_count);
  for (int i = 0; i < tree->tip_count; ++i) {
    ENTRY entry;
    entry.key = tree->nodes[i]->label;
    entry.data = &tree->nodes[i]->clv_index;
    hsearch(entry, ENTER);
  }
  return tree;
}

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

pll_partition_t *createPartition(const pll_msa_t *msa)
{

  pll_partition_t *partition = pll_partition_create(msa->count,
      msa->count - 1, // inner node count
      4, // states
      (unsigned int)(msa->length),
      1,
      2 * msa->count - 2, // branch count
      1,
      msa->count - 1, // inner node count
      PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_SITE_REPEATS);
  if (!partition) {
    cerr << "[Error] Could not create a partition" << endl;
  }
  return partition;
}

void createPartitions(const vector<pll_msa_t *> &msas,
  vector<pll_partition_t *> &partitions)
{
  for (auto msa: msas) {
    partitions.push_back(createPartition(msa));
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

  pll_utree_t * tree = loadTree(newickFile);
  pll_msa_t *fullMSA = loadMSA(msaFile);
  vector<pll_msa_t *> partitionnedMSAs;
  partitionMSA(fullMSA, partFile, partitionnedMSAs);  
  vector<pll_partition_t *> partitions;
  createPartitions(partitionnedMSAs, partitions); 
  cout << "Successfuly wrote output into " << outputFile << endl;
  return 0;
}
