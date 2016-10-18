for i in `ls *mesh.fds`
do
	echo $i
	sed -f change.cmd $i > tmp
	mv tmp $i
done
