#ifndef _REPEATS_COUNTER_DATASET_HPP_
#define _REPEATS_COUNTER_DATASET_HPP_

#include <string> 

namespace RepeatsCounter {

class Dataset {
public:
  Dataset(const std::string &fastaFile,
      const std::string &newickFile);

  int getRepeatsClassesNumber() const;
};



} // namespace RepeatsCounter

#endif
