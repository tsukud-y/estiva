default: stokes

stokes: ;
	(cd plot; cc -c *.c -I../)
	cc *.c plot/*.o -L/usr/local/lib/ -lm -lf2c -o stokes -lestiva

clean: ;
	/bin/rm -f *.o plot/*.o stokes
