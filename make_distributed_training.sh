rm CMakeCache.txt
cmake .
make clean
make VERBOSE=1 distributed_training_main -j
make VERBOSE=1 distributed_training_moderator -j 
