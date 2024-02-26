python3.10 01
grep "Atomic Numbers:" output.txt > Atno_N.txt
sed -i "s+Atomic Numbers: tensor++g" Atno_N.txt
sed -i "s+(++g" Atno_N.txt
sed -i "s+)++g" Atno_N.txt
sed -i '/Atomic Numbers/d' output.txt
sed -i "s+Representation: ++g" output.txt
python3.10 02 > desc_N.dat
rm output.txt Atno_N.txt
python3 03
