SHELL = /bin/sh
export ROBOT = /usr/bin/perl -I/home/bp/oofem-doc/www/en /home/bp/oofem-doc/www/robot.pl

html: 
	cd elementlibmanual; make html
	cd oofemInput; make html
	cd paper-actap; make html
	cd extractorInput; make html
	cd matlibmanual; make html
	cd programmer; make html
	cd refman; make html

ps:
	cd elementlibmanual; make ps
	cd oofemInput; make ps
	cd paper-actap; make ps
	cd extractorInput; make ps
	cd matlibmanual; make ps
	cd programmer; make ps

pdf:
	cd elementlibmanual; make pdf
	cd oofemInput; make pdf
	cd paper-actap; make pdf
	cd extractorInput; make pdf
	cd matlibmanual; make pdf
	cd programmer; make pdf

www-install:
ifdef INSTALLDIR
	export INSTALLDIR
	cd elementlibmanual; make www-install
	cd oofemInput; make www-install
	cd paper-actap; make www-install
	cd extractorInput; make www-install
	cd matlibmanual; make www-install
	cd programmer; make www-install
	cd refman; make www-install
else
	echo "INSTALLDIR not defined"
endif


clean:
	cd elementlibmanual; make clean
	cd oofemInput; make clean
	cd paper-actap; make clean
	cd extractorInput; make clean	
	cd matlibmanual; make clean
	cd programmer; make clean
	cd refman; make clean
