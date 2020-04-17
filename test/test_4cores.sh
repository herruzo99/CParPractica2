# with set -e if any diff fails the scripts ends
set -e
mkdir -p output_modificado
# 150 Standard output_modificado with different arguments
for i in {1..20}
do
echo "Fase $(($i))"
mpiexec -n $(($i*13%3 + 1)) ../evolution 200 25 1000 1000 1 10 $((4424+$i*127)) 6534 2461 $(($i*19059 % 4000 +1)) 5 5 10 10 10 10 | sed --expression='$!d' >output_modificado/Iter$(($i))Part1Test; diff  output_original/Iter$(($i))Part1Test output_modificado/Iter$(($i))Part1Test
mpiexec -n $(($i*4%3 + 1)) ../evolution 30 30 1500 1000 1 7 4624 $((9334+$i*7)) 961 10 | sed --expression='$!d' >output_modificado/Iter$(($i))Part2Test; diff  output_original/Iter$(($i))Part2Test output_modificado/Iter$(($i))Part2Test
mpiexec -n $(($i*5%3 + 1)) ../evolution 40 10 200 1000 1 2 4468 57764 $((2461+$i*13)) $(($i*4741 % 40 +1)) | sed --expression='$!d' >output_modificado/Iter$(($i))Part3Test; diff  output_original/Iter$(($i))Part3Test output_modificado/Iter$(($i))Part3Test
mpiexec -n $(($i*6%3 + 1)) ../evolution 200 10 100 1000 2 4 $((44+$i*37)) 55534 $((9361+$i*56)) 100 1 1 50 1 3 1 | sed --expression='$!d' >output_modificado/Iter$(($i))Part4Test; diff  output_original/Iter$(($i))Part4Test output_modificado/Iter$(($i))Part4Test
mpiexec -n $(($i*11%3 + 1)) ../evolution 5 500 75 1000 1 11 $((4324+$i*8)) $((576+$i*3)) $((461+$i*11)) $(($i*742 % 80 +1)) | sed --expression='$!d' >output_modificado/Iter$(($i))Part5Test; diff  output_original/Iter$(($i))Part5Test output_modificado/Iter$(($i))Part5Test
done
# Special Cases
#  A lot of iterations
echo "Fase 1x1"
mpiexec -n 3 ../evolution 1 1 1000000 1000 5 10 444324 5776534 9542462 5 | sed --expression='$!d' >output_modificado/OneTileTest; diff -q output_original/OneTileTest output_modificado/OneTileTest
#  A lot of cells on a small space
echo "Fase 2x2"
mpiexec -n 2 ../evolution 2 2 500000 1000 10 10 434324 5776533 9542462 5 | sed --expression='$!d' >output_modificado/2x2TileTest; diff -q output_original/2x2TileTest output_modificado/2x2TileTest
#  A lot of cells on a big space
echo "Fase Grande"
mpiexec -n 4 ../evolution 1000 1000 100 1000 5 10 444324 5776534 9542462 500 | sed --expression='$!d' >output_modificado/BigTest; diff -q output_original/BigTest output_modificado/BigTest
#  A lot of cells without food
echo "Fase Muertos"
mpiexec -n 1 ../evolution 1000 1000 100 1000 0 0 444324 5776534 9542462 2000000 | sed --expression='$!d' >output_modificado/AllDeadTest;  diff -q output_original/BigTest output_modificado/BigTest
#  A lot of food in one spot
echo "Fase Mucha comida en un punto"
mpiexec -n 2 ../evolution 1000 1000 500 100000 0 0 444324 5776534 9542462 2000  500 500 40 40 10 3 | sed --expression='$!d' >output_modificado/SpotFoodTest;  diff -q output_original/SpotFoodTest output_modificado/SpotFoodTest

echo "Todos los test completados con éxito"
