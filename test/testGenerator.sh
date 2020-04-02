#Generate original executable
make -C.. clean
make -C../original
mv ../original/evolution ../
#TestGenerator es un script para genear las pruebas
mkdir -p testBase
for i in {1..30}
do
../evolution 20 20 1000 1000 1 1 444324$i 25776534 9542461 1 | sed --expression='$!d' >testBase/Seed1Test$i
../evolution 20 20 1000 1000 1 1 444324 25776534$i 9542461 1 | sed --expression='$!d' >testBase/Seed2Test$i
../evolution 20 20 1000 1000 1 1 444324 25776534 9542461$i 1 | sed --expression='$!d' >testBase/Seed3Test$i
../evolution 20 20 1000 1000 1 1 444324$i 25776534 9542461$i 1 | sed --expression='$!d' >testBase/Seed13Test$i
../evolution 20 20 1000 1000 1 1 444324$i 25776534$i 9542461$i 1 | sed --expression='$!d' >testBase/Seed123Test$i
done
../evolution 1 1 1000000 1000 5 10 444324 25776534 9542462 5 | sed --expression='$!d' >testBase/OneTileTest
../evolution 2 2 500000 1000 10 10 434324 25776533 9542462 5 | sed --expression='$!d' >testBase/2x2TileTest
../evolution 1000 1000 100 1000 5 10 444324 25776534 9542462 5 | sed --expression='$!d' >testBase/BigTest
../evolution 1000 1000 100 1000 0 0 444324 25776534 9542462 2000000 | sed --expression='$!d' >testBase/AllDeadTest
#Restore executable
make -C.. clean
make -C..