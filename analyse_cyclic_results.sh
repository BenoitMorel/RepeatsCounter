mkdir -p build
cd build
cmake ..
make
cd ..
build/RepeatsCounter datasets/59.repeats datasets/59.cyclic.txt
echo ""
echo ""
build/RepeatsCounter datasets/128.repeats datasets/128.cyclic.txt
echo ""
echo ""
build/RepeatsCounter datasets/404.repeats datasets/404.cyclic.txt

