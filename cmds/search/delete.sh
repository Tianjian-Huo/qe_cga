for dir in `ls -d */|grep '/$'`
do
	cd $dir
		rm dmol.tpdensk
		rm dmol.tpotl
		echo "$dir finished"
	cd ..
done
