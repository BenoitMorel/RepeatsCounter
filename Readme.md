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
benoit:~/github/RepeatsCounter/build$ ./RepeatsCounter ../example/small_dataset/small.repeats ../example/small_dataset/small.distribution
Parsing repeats file ../example/small_dataset/small.repeats...
Parsing distribution file ../example/small_dataset/small.distribution...
*******************************************
**  Analyse with one single core   **
*******************************************
UniqueCore RCC: 26   sites: 10   partitions: 2
Worst RCC:    26
Worst RCC * cores:  26
Sum of RCCs:    26
*******************************************
**  Analyse with multiple cores  **
*******************************************
CPU1 RCC: 13   sites: 5  partitions: 1
CPU2 RCC: 14   sites: 5  partitions: 2
Worst RCC:    14
Worst RCC * cores:  28
Sum of RCCs:    27
```

The first analysis is the sequential RCC (with only one core, not parallelization, no split). You can see it like the ultimate lower bound, that you can not reach.
The second analysis is run on your distribution file.
- The worst RCC is the target to minimize, for a given repeats file and a given number of cores. 
- The sum of RCCs is the sum over the cores of there respective RCC. The less you split repeats, the closer it will be to the unique-core RCC.
- The Worst RCC * cores will get closer to the sum of RCC if your load balance gets better.
To acheive this, you need to reduce the total cost while keeping the load balance among cores as good as possible.
This example is useful to undertand files syntaxes, but you should better work on the files in datasets.tar.gz.



## Generating repeats file (for teachers)

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


