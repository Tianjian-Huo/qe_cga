#!/bin/bash
for i in `ls -d */|grep '/$'`:
do
	cd $i
		mv dmol.car dmol_car_pro
                rm *.car
                mv dmol_car_pro dmol.car
                
	cd ..
done
