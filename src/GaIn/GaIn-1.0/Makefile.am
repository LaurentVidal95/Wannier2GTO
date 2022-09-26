EXTRA_DIST = doc

SUBDIRS = \
	src \
	tests

DISTCLEANFILES = \
	lib/libGaIn.a \
	include/GaIn.FC.h \
	include/GaIn.h 
		
all-am: 
	mkdir -p lib; mv -f src/libGaIn.a lib/.;
	mkdir -p include; cp -f $(top_srcdir)/src/GaIn.h include/.
	grep "FC_FUNC(" src/GaIn.config.h > include/GaIn.FC.h
	sed -i "s/FC_FUNC/GAIN_FC_FUNC/g" include/GaIn.FC.h

