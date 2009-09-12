default: stokes

stokes: ;
	(cd solver; cc -c *.c -I../ -I/usr/local/include/)
	(cd plot; cc -c *.c -I../)
	cc *.c solver/*.o plot/*.o -L/usr/local/lib/ -lm -lf2c -o stokes -lestiva

clean: ;
	/bin/rm -f *.o solver/*.o plot/*.o stokes
