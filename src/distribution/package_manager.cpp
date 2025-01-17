#include "distribution/package_manager.hpp"

PackageManager* PackageManager::manager = nullptr;
std::mutex PackageManager::mtx;