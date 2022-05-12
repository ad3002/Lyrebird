# Lyrebird

Raw reads based satellite DNA extraction pipeline

### Prerequisites

1. Currently only Linux supported and requeired some skills in software installation and running. I'm going to pack it into Docker image in the near future then it will be easy to use. And you can always contact me with ad3002@gmail.com and I will run your data.

2. Build aindex:

```bash
cd aindex
cd external
git clone https://github.com/ot/emphf.git
cd emphf
cmake .
make
mkdir bin
g++ Compute_index.cpp  -std=c++11 -pthread -O3 debrujin.hpp debrujin.cpp read.hpp read.cpp kmers.hpp kmers.cpp settings.cpp settings.hpp hash.hpp hash.cpp && mv a.out bin/compute_index.exe
g++ -c -std=c++11 -fPIC python_wrapper.cpp -o python_wrapper.o && g++ -c -std=c++11 -fPIC kmers.cpp kmers.hpp debrujin.cpp debrujin.hpp hash.cpp hash.hpp read.cpp read.hpp settings.hpp settings.cpp && g++ -shared -Wl,-soname,python_wrapper.so -o python_wrapper.so python_wrapper.o kmers.o debrujin.o hash.o read.o settings.o
cd ../..
```

3. You need MongoDB (see https://www.mongodb.com/docs/manual/administration/install-on-linux/) and python driver for it:

```
pip3 install pymongo
```

4. Other python packages:

```
pip3 install simplejson
```

4. You need jellyfish or any other kmer counter, you need tab-delimeted file with kmer and kmer_frequency. Jellyfish can be obtained here https://github.com/gmarcais/Jellyfish/releases 

### Usage

1. IMPORTANT! Carefully remove all adapters sequences from raw reads with any tool of your choice. Otherwise they will produce multiple false positive predictions. 

2. Compute tab-delimted file with Jellyfish (for details see JF docs):

```bash
jellyfish count -m 23 -t THREAD -s MEMORY -C -o PREFIX.23.jf2 INPUT_FILES
jellyfish dump -t -c -o PREFIX.23.dat PREFIX.23.jf2
```

3. Compute index with dat file:

```bash
sort -k2nr PREFIX.23.dat > PREFIX.23.sdat
cut -f1 PREFIX.23.dat > PREFIX.23.kmers"
./aindex/emphf/compute_mphf_seq.exe PREFIX.23.kmers PREFIX.23.pf"
./aindex/bin/compute_index.exe PREFIX.23.dat PREFIX.23.pf PREFIX.23 THREADS 0
```

4. Run Lyrebird:

Lyrebird searches tandem repeats circles in de Brujin graph, so you need precomputed and sorted file with kmer ferquencies:

```
lyrebird.py -a PREFIX.23 -i PREFIX.23.sdat -o OUTPUT_PREFIX -p SAT_PREFIX --coverage 1

```

SAT_PREFIX - desired prefix for satellite DNA name