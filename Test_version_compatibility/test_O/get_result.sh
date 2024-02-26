python3.10 01
grep "Atomic Numbers:" output.txt > Atno_O.txt
sed -i "s+Atomic Numbers: tensor++g" Atno_O.txt
sed -i "s+(++g" Atno_O.txt
sed -i "s+)++g" Atno_O.txt
sed -i '/Atomic Numbers/d' output.txt
sed -i "s+Representation: ++g" output.txt
python3.10 02 > desc_O.dat
rm output.txt Atno_O.txt
python3 03
