.SUFFIXES:  .d

CSRC = $(shell ls *.cc)
INC =  $(shell ls *.h )

OBJ = $(CSRC:.cc=.o)

all:: $(OBJ)

.PHONY: clean

# this includes the dependecy files
-include $(OBJ:.o=.d)

# The next line gets the dependences
.cc.d:
	set -e; $(CPP) -MM $(CPP_FLAGS) $< \
	| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
	[ -s $@ ] || rm -f $@

.cc.o:
	${CPP} ${CPP_FLAGS}  -c $(CPP_FLAGS)  $<


clean:
	rm *~  *.o *.d
