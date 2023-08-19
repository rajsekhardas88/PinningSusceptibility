./logIndices.sh 200 1.4 5000000 | awk '{if($1>0)print $0}' > tmp1.txt
tail -1 tmp1.txt > a.dat
awk '{print $1+1}' a.dat > b.dat
cat tmp1.txt b.dat > tmp.txt
rm a.dat b.dat tmp1.txt
