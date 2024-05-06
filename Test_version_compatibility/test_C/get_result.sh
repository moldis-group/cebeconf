python3.10 01
grep "Atomic Numbers:" output.txt > Atno_C.txt
sed -i "s+Atomic Numbers: tensor++g" Atno_C.txt
sed -i "s+(++g" Atno_C.txt
sed -i "s+)++g" Atno_C.txt
sed -i '/Atomic Numbers/d' output.txt
sed -i "s+Representation: ++g" output.txt
python3.10 02 > desc_C.dat
rm output.txt Atno_C.txt
python3 03
