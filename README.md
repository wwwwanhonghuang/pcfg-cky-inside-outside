!!! Current Result are incorrect, fixing.
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

### Run
`$ ./main <pcfg-grammar-file-path> <input-sequence-file-path>`

## Reference
[1] Itiroo Sakai, “Syntax in universal translation”. In Proceedings 1961 International Conference on Machine Translation of Languages and Applied Language Analysis, Her Majesty’s Stationery Office, London, p. 593-608, 1962.

[2] Grune, Dick (2008). Parsing techniques : a practical guide (2nd ed.). New York: Springer. p. 579. ISBN 978-0-387-20248-8.

[3] Stratos, Karl. “The Inside-Outside Algorithm.” (2012).
