


### primarily is used to move the pickle files to main_pickle/ directory

mkdir main_pickle
for i in {0..5}; do
	mv merge_out/$i-out/main.pickle merge_out/$i-out/main-$i.pickle
	cp merge_out/$i-out/main-$i.pickle main_pickle
	mv merge_out/$i-out/main-$i.pickle merge_out/$i-out/main.pickle
done
