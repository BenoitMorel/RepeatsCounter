# RepeatsCounter

RepeatsCounter is a tool that computes the RCC (repeats class count) from a repeats file and a data distribution file.

It will be used to evaluate the solutions of the students doing the Bioinformatics practical of 2018.

If you have any question or if you find some bug, please write me a mail at benoit.morel@h-its.org

## Building RepeatsCounter

You will need cmake to compile the project.

Run the following:
```
cd RepeatsCounter
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


# RDDA
A quick implementation of RDDA is provided 

Compilation:
```
cd RDDA
mkdir build
cd build
cmake ..
make
```

Running:
```
./rdda repeats_file cores output_file
```


## Generating repeats file from MSA and partition file (for teachers)

@Students: you don't need to read this paragraph :-)

To generate a repeats file from a tree and a partionned MSA, you need to clone the project with --recursive option (to download libpll as well)
Then:
```
cd tools/MSAConverter
mkdir build
cd build
cmake ..
make
```

and

```
./convert msa_file partition_file newick_file output_repeats_file
```


