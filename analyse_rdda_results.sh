cd RepeatsCounter
mkdir -p build
cd build
cmake ..
make
cd ../..
cores=10
RepeatsCounter/build/RepeatsCounter datasets/59.repeats results/59_${cores}cores.rdda.txt
echo ""
echo ""
RepeatsCounter/build/RepeatsCounter datasets/128.repeats results/128_${cores}cores.rdda.txt
echo ""
echo ""
RepeatsCounter/build/RepeatsCounter datasets/404.repeats results/404_${cores}cores.rdda.txt

cores=2
RepeatsCounter/build/RepeatsCounter datasets/59.repeats results/59_${cores}cores.rdda.txt
echo ""
echo ""
RepeatsCounter/build/RepeatsCounter datasets/128.repeats results/128_${cores}cores.rdda.txt
echo ""
echo ""
RepeatsCounter/build/RepeatsCounter datasets/404.repeats results/404_${cores}cores.rdda.txt

cd ../..
