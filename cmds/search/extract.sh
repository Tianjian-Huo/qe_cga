#!/bin/bash
n = 0
for i in `ls -d */|grep '/$'`:
do
    cd $i
        cp ../extract_car_opted.py .
        tail -28 dmol.outmol|head -16 > str.txt
        first_str=`head -1 str.txt | cut -c  13`
        if [ $first_str == 1 ]
        then
            n=$(($n+1))
            python3 extract_car_opted.py
            mv test.car $n.car
            rm str.txt extract_car_opted.py
            echo "        "
            echo " ------------------ $n finished ------------------"
            echo "        "
        else
            echo "optmization failed"
        fi
    cd ..
done
