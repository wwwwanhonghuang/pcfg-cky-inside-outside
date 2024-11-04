#ifndef H_DATASET_HELPER
#define H_DATASET_HELPER

#include <vector>
#include <cstddef>
#include <stdexcept>
#include <cstdint>
#include <random>
#include <algorithm>

void split_dataset(const std::vector<std::vector<uint32_t>>& sentences,
                  std::vector<std::vector<uint32_t>>& train_set,
                  std::vector<std::vector<uint32_t>>& valid_set,
                  double train_fraction);

#endif
