#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <search.h>
#include "pll.h"
#include <memory>

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

pll_msa_t *loadMSAFasta(const string &fastaFile)
{
  auto reader = pll_fasta_open(fastaFile.c_str(), pll_map_fasta);
  if (!reader) {
    cerr << "[ERROR] Cannot open " << fastaFile << endl;
    return 0;
  }
  char * head;
  long head_len;
  char *seq;
  long seq_len;
  long seqno;
  int length;
  vector<pair<string, string> > sequences;
  while (pll_fasta_getnext(reader, &head, &head_len, &seq, &seq_len, &seqno)) {
    string sequence(seq, seq_len);
    string label(head, head_len);
    sequences.push_back(pair<string, string>(sequence, label));
  }
  pll_msa_t * msa = (pll_msa_t *)malloc(sizeof(pll_msa_t));
  msa->count = sequences.size();
  msa->length = sequences[0].first.size();
  msa->label = (char**)malloc(msa->count * sizeof(char *));
  msa->sequence = (char**)malloc(msa->count * sizeof(char *));
  for (unsigned int i = 0; i < msa->count; ++i) {
    auto &sequence = sequences[i].first;
    auto &label = sequences[i].second;
    msa->sequence[i] = (char *)malloc(sizeof(char) * (sequence.size() + 1));
    msa->label[i] = (char *)malloc(sizeof(char) * (label.size() + 1));
    memcpy(msa->sequence[i], sequence.c_str(), sequence.size());
    memcpy(msa->label[i], label.c_str(), label.size());
    msa->sequence[i][sequence.size()] = 0;
    msa->label[i][label.size()] = 0;
  }
  /*
  weights = pll_compress_site_patterns(buffer, map, count, &length);
  if (!weights) 
    throw LibpllException("Error while parsing fasta: cannot compress sites");
  for (unsigned int i = 0; i < count; ++i) {
    sequences[i]->len = length;
  }
  */
  pll_fasta_close(reader);
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
      memcpy(submsa->sequence[i], msa->sequence[i] + first, submsa->length);
    }
    partitionnedMSAs.push_back(submsa);
  }
}

pll_partition_t *createPartition(const pll_msa_t *msa)
{
  const int RATES = 1;
  pll_partition_t *partition = pll_partition_create(msa->count,
      msa->count - 1, // inner node count
      4, // states
      (unsigned int)(msa->length),
      1,
      2 * msa->count - 2, // branch count
      RATES,
      msa->count - 1, // inner node count
      PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_SITE_REPEATS);
  if (!partition) {
    cerr << "[Error] Could not create a partition" << endl;
  }
  
  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };
  double subst_params[6] = {1,1,1,1,1,1};
  double rate_cats[RATES] = {1.0};
  pll_set_frequencies(partition, 0, frequencies);
  pll_set_subst_params(partition, 0, subst_params);
  pll_set_category_rates(partition, rate_cats);

  for (int i = 0; i < msa->count; ++i)
  {
    ENTRY query;
    query.key = msa->label[i];
    ENTRY *found = hsearch(query,FIND);
    if (!found) {
      cerr << "[Error] Cannot find node " << msa->label[i] << " in the newick tree " << endl;
    }
    unsigned int tipIndex = *((unsigned int *)(found->data));
    pll_set_tip_states(partition, tipIndex, pll_map_nt, msa->sequence[i]);
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

static int cb_full_traversal(pll_unode_t * node)
{
  return 1;
}

void updatePartials(vector<pll_partition_t *> &partitions, pll_utree_t *utree)
{
  for (auto partition: partitions) {
    pll_unode_t **travbuffer = (pll_unode_t **)malloc(partition->nodes * sizeof(pll_unode_t *));
    unsigned int traversalSize;
    pll_unode_t *root = utree->nodes[0];
    pll_utree_traverse(root,
          PLL_TREE_TRAVERSE_POSTORDER,
          cb_full_traversal,
          travbuffer,
          &traversalSize);
    int branches = 2 * utree->tip_count - 2;
    double *branch_lengths = (double *)malloc(branches * sizeof(double));
    unsigned int *matrix_indices = (unsigned int *)malloc(branches * sizeof(unsigned int));
    pll_operation_t *operations = (pll_operation_t *)malloc(utree->inner_count *
          sizeof(pll_operation_t));
    unsigned int matrix_count, ops_count;
    pll_utree_create_operations(travbuffer,
        traversalSize,
        branch_lengths,
        matrix_indices,
        operations,
        &matrix_count,
        &ops_count);
    for (int i = 0; i < ops_count; ++i) {
      pll_update_repeats(partition, operations + i);
    }
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
  pll_msa_t *fullMSA = loadMSAFasta(msaFile);
  vector<pll_msa_t *> partitionnedMSAs;
  partitionMSA(fullMSA, partFile, partitionnedMSAs);  
  vector<pll_partition_t *> partitions;
  createPartitions(partitionnedMSAs, partitions); 
  updatePartials(partitions, tree); 
  cout << "Successfuly wrote output into " << outputFile << endl;
  return 0;
}
