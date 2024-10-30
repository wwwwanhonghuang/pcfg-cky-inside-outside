#include "dataset_helper.hpp"
void split_dataset(const std::vector<std::vector<uint32_t>>& sentences,
                  std::vector<std::vector<uint32_t>>& train_set,
                  std::vector<std::vector<uint32_t>>& valid_set,
                  double train_fraction) {
    if (train_fraction <= 0.0 || train_fraction >= 1.0) {
        throw std::invalid_argument("train_fraction must be between 0 and 1.");
    }

    std::vector<std::vector<uint32_t>> shuffled_sentences = sentences;

    std::random_device rd;
    std::mt19937 generator(rd());
    std::shuffle(shuffled_sentences.begin(), shuffled_sentences.end(), generator);

    int total_sentences = shuffled_sentences.size();
    int split_point = static_cast<size_t>(total_sentences * train_fraction);

    train_set.assign(shuffled_sentences.begin(), shuffled_sentences.begin() + split_point);
    valid_set.assign(shuffled_sentences.begin() + split_point, shuffled_sentences.end());
}