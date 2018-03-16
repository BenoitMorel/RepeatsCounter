#include <iostream>
#include <string>

#include "pll.h"

using namespace std;

void loadFullMSA() {
  pll_phylip_t * fd = pll_phylip_open("", 0);
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

  cout << "Not implemented" << endl; 

  return 0;
}
