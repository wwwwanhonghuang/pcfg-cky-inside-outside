!!! Current Result may be incorrect. The appliaction is under testing.
# CKY Parser and Inside-outside Algorithm for PCFG
This repository contains a C++ implementation of CYK parser[1, 2] and inside-outside algorithm[3] for
estimating a PCFG.
Currently, CPU-only version is under testing.
CUDA version is partially implemented.

Current implementation may contain errors.
## Usages
### Prepare the PCFG Definition File
1. A PCFG definiton file contain multiple line, with each line define a PCFG rule. 
2. A PCFG rule in the PCFG definiton file has the form of `A->B C [possibility]` or `A->B [possibility]`
3. In a rule's right side, currently support $2$ symbols at most.
4. The rule is not constrained to be a CNF form. i.e., **for $B$ or $C$, it can either be a non-termiante or a terminate**.
5. A termiante should quote by '\'\''. Example `A->'w_1' B`, `A->B 'w_1'`, `A->'w_1' 'w_2'`, where `w_1` and `w_2` are termiantes.  
6. Example definition files are provided in the root of this repository, with suffix of `.pcfg`

### Prepare the Input Sequence File
1. An input sequence file may contains multiple lines.
2. Each line is a word sequence. Each word separated by a space.
3. Currently, each line's works will be concated into one sequence. The parser will analyze a sequence of [line1_words  line2_words line3_words ...].

### Make
`$ make` 

### Configuration
You may need edit the `config.yaml` file at the root of this repository, to setup proper PCFG estimation 
configurations.

Example
``` yaml
main:
  grammar_file: "data/grammar.pcfg"  # Path to the grammar file used for parsing sentences.
  input: "data/sentences_converted.txt"  # Path to the file with sentences to be processed.
  log_intervals: 20000  # Interval for logging rule probabilities during each epoch. Ignored if batch updates are disabled.
  log_path: "data/logs/"  # Directory where log files will be saved.
  
  log_f:
    enabled: true  # Toggle to enable logging of the 'f' variable in the inside-outside algorithm.
    intervals: 20000  # Interval for logging the 'f' variable during training.

  log_warning_in_training: false  # Enable/disable logging of warnings during training.

  batch_size_for_parameter_update: -1  # If positive, updates grammar rule probabilities at this interval during each epoch. If set to -1, updates are performed only at the end of each epoch. Midway updates can impact convergence.

  split_data:
    enabled: true  # Enable/disable splitting data into training and validation sets.
    val_dataset_path: "data/validate_sentences.txt"  # Path to the validation dataset.
    train_dataset_path: "data/train_sentences.txt"  # Path to the training dataset.
    train_fraction: 0.8  # Fraction of data used for training; the rest is used for validation.

  n_epochs: 10  # Total number of epochs for training.
  limit_n_sentences: -1  # Limits the number of sentences processed from the input file. If set to -1, all sentences are used.
```
### Run
`$ ./build/bin/main_executable`

## Reference
[1] Itiroo Sakai, “Syntax in universal translation”. In Proceedings 1961 International Conference on Machine Translation of Languages and Applied Language Analysis, Her Majesty’s Stationery Office, London, p. 593-608, 1962.

[2] Grune, Dick (2008). Parsing techniques : a practical guide (2nd ed.). New York: Springer. p. 579. ISBN 978-0-387-20248-8.

[3] Stratos, Karl. “The Inside-Outside Algorithm.” (2012).
