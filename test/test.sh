set -e
mkdir -p tests
# 150 Standard tests with diff -qerent seeds
for i in {1..30}
do
../evolution 20 20 1000 1000 1 1 444324$i 25776534 9542461 1 | sed --expression='$!d' >tests/Seed1Test$i; diff -q testBase/Seed1Test$i tests/Seed1Test$i
../evolution 20 20 1000 1000 1 1 444324 25776534$i 9542461 1 | sed --expression='$!d' >tests/Seed2Test$i; diff -q testBase/Seed1Test$i tests/Seed1Test$i
../evolution 20 20 1000 1000 1 1 444324 25776534 9542461$i 1 | sed --expression='$!d' >tests/Seed3Test$i; diff -q testBase/Seed1Test$i tests/Seed1Test$i
../evolution 20 20 1000 1000 1 1 444324$i 25776534 9542461$i 1 | sed --expression='$!d' >tests/Seed13Test$i; diff -q testBase/Seed1Test$i tests/Seed1Test$i
../evolution 20 20 1000 1000 1 1 444324$i 25776534$i 9542461$i 1 | sed --expression='$!d' >tests/Seed123Test$i; diff -q testBase/Seed1Test$i tests/Seed1Test$i
done
# Special Cases
#  A lot of iterations
../evolution 1 1 1000000 1000 5 10 444324 25776534 9542462 5 | sed --expression='$!d' >tests/OneTileTest; diff -q testBase/OneTileTest tests/OneTileTest
#  A lot of cells on a small space
../evolution 2 2 500000 1000 10 10 434324 25776533 9542462 5 | sed --expression='$!d' >tests/2x2TileTest; diff -q testBase/2x2TileTest tests/2x2TileTest
#  A lot of cells on a big space
../evolution 1000 1000 100 1000 5 10 444324 25776534 9542462 5 | sed --expression='$!d' >tests/BigTest; diff -q testBase/BigTest tests/BigTest
#  A lot of cells without food
../evolution 1000 1000 100 1000 0 0 444324 25776534 9542462 2000000 | sed --expression='$!d' >testBase/AllDeadTest;  diff -q testBase/BigTest tests/BigTest


echo "Todos los test completados con éxito"
