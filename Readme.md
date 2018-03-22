# RepeatsCounter

RepeatsCounter is a tool that computes the RCC (repeats class count) from a repeats file and a data distribution file.

It will be used to evaluate the solutions of the students doing the Bioinformatics practical of 2018.

If you have any question or if you find some bug, please write me a mail at benoit.morel@h-its.org

## Building RepeatsCounter

You will need cmake to compile the project.

Run the following:
```
mkdir build
cd build
cmake ..
make
```
## Running the program

From the build directory:
```
./RepeatsCounter repeats_file distribution_file
```

For instance, from the build directory:
```
./RepeatsCounter ../example/small_dataset/small.repeats ../example/small_dataset/small.distribution 
Parsing repeats file ../example/small_dataset/small.repeats...
Parsing distribution file ../example/small_dataset/small.distribution...
Core CPU1 cost:   13
Core CPU2 cost:   14
Worst core cost:  14
Total cost:     27
```

Cost means RCC.
The worst core cost is the target to minimize, for a given repeats file and a given number of cores. 
To acheive this, you need to reduce the total cost while keeping the load balance among cores as good as possible.
This example is useful to undertand files syntaxes, but you should better work on the files in datasets.tar.gz.

## Generating repeats file (for teachers)

@Students: you don't need to read this paragraph :-)

To generate a repeats file from a tree and a partionned MSA, you need to clone the project with --recursive option (to download libpll as well)
Then:
```
cd MSAConverter
mkdir build
cd build
cmake ..
make
```

and

```
./convert msa_file partition_file newick_file output_repeats_file
```


