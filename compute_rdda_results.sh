cd RDDA/
mkdir -p build
cd build
cmake ..
make
cd ../../
mkdir -p results
cores=10
RDDA/build/rdda datasets/59.repeats ${cores} results/59_${cores}cores.rdda.txt
RDDA/build/rdda datasets/128.repeats ${cores} results/128_${cores}cores.rdda.txt
RDDA/build/rdda datasets/404.repeats ${cores} results/404_${cores}cores.rdda.txt

