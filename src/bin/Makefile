default : all

all : ;
	(cd xyc2msh;make)
	(cd xmsh;make)
	(cd xmesh;make)
	(cd unestiva;make)

clean : ;
	(cd xyc2msh; make clean)
	(cd xmsh; make clean)
	(cd xmesh; make clean)
	(cd unestiva; make clean)

install : ;
	cp xyc2msh/xyc2msh xmsh/xmsh xmesh/xmesh unestiva/unestiva /usr/local/bin
