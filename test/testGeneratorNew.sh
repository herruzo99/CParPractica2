#Generate original executable
make -C.. clean
make -C../original
mv ../original/evolution ../
#TestGenerator es un script para genear las pruebas
mkdir -p testBaseNew
for i in {1..20}
do
echo "Fase $(($i))"
../evolution 200 25 1000 1000 1 10 $((4424+$i*127)) 6534 2461 $(($i*19059 % 4000 +1)) | sed --expression='$!d' >testBaseNew/Seed1Test$i
../evolution 30 30 1500 1000 1 7 4624 $((9334+$i*7)) 961 10 | sed --expression='$!d' >testBaseNew/Seed2Test$i
../evolution 40 10 200 1000 1 2 4468 57764 $((2461+$i*13)) $(($i*4741 % 40 +1)) | sed --expression='$!d' >testBaseNew/Seed3Test$i
../evolution 200 10 100 1000 2 4 $((44+$i*37)) 55534 $((9361+$i*56)) 100 | sed --expression='$!d' >testBaseNew/Seed13Test$i
../evolution 5 500 75 1000 1 11 $((4324+$i*8)) $((576+$i*3)) $((461+$i*11)) $(($i*742 % 80 +1)) | sed --expression='$!d' >testBaseNew/Seed123Test$i
done
echo "Fase 1x1"
../evolution 1 1 1000000 1000 5 10 444324 5776534 9542462 5 | sed --expression='$!d' >testBaseNew/OneTileTest
echo "Fase 2x2"
../evolution 2 2 500000 1000 10 10 434324 5776533 9542462 5 | sed --expression='$!d' >testBaseNew/2x2TileTest
echo "Fase Grande"
../evolution 1000 1000 100 1000 5 10 444324 5776534 9542462 500 | sed --expression='$!d' >testBaseNew/BigTest
echo "Fase Muertos"
../evolution 1000 1000 100 1000 0 0 444324 5776534 9542462 2000000 | sed --expression='$!d' >testBaseNew/AllDeadTest
#Restore executable
make -C.. clean
make -C..
