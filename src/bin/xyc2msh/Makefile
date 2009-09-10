default : xyc2msh

xyc2msh : ;
	(cd oldestiva;make)
	gcc *.c -Ioldestiva oldestiva/*.o -o xyc2msh

clean : ;
	rm -f *.o xyc2msh
	(cd oldestiva; make clean)
