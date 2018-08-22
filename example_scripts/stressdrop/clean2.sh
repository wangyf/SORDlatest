
#cp ../../../static_stressdrop/matlabRUPGENv1/slip.bin bin/
#python strdrop.py -f
cp bin/slip.bin strdrop2/in
rm strdrop2/out/* strdrop2/stats/* strdrop2/debug/* strdrop2/prof/*
cp src/sord-sO strdrop2/
sh strdrop2/run.sh