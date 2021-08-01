# LIPP

Thanks for all of you for your interest in our work.

LIPP, the Updatable Learned Index with Precise Positions, is a brand new learned index to provide efficient support for a set of index operations. 

This project contains the code of LIPP and welcomes contributions or suggestions.

## Compile & Run

```bash
mkdir build
cd build
cmake ..
make
```

Run example:

```bash
./build/example
./build/example_without_fmcd
./build/example_bulk_load
```

## Usage

See `src/examples/*` and `src/core/lipp.h` for detail.

## Publication

> Wu J, Zhang Y, Chen S, et al. Updatable Learned Index with Precise Positions. PVLDB, 14(8): 1276-1288 (2021)

```
@article{DBLP:journals/pvldb/WuZCCWX21,
  author    = {Jiacheng Wu and
               Yong Zhang and
               Shimin Chen and
               Yu Chen and
               Jin Wang and
               Chunxiao Xing},
  title     = {Updatable Learned Index with Precise Positions},
  journal   = {Proc. {VLDB} Endow.},
  volume    = {14},
  number    = {8},
  pages     = {1276--1288},
  year      = {2021},
  url       = {http://www.vldb.org/pvldb/vol14/p1276-wu.pdf},
  timestamp = {Wed, 14 Jul 2021 12:34:06 +0200},
  biburl    = {https://dblp.org/rec/journals/pvldb/WuZCCWX21.bib},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```
