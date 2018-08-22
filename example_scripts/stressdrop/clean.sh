
#cp ../../../static_stressdrop/matlabRUPGENv1/slip.bin bin/
#python strdrop.py -f
cp bin/slip.bin strdrop/in
rm strdrop/out/* strdrop/stats/* strdrop/debug/* strdrop/prof/*
cp src/sord-sO strdrop/
sh strdrop/run.sh