TARGET = crabber
MODULES = $(patsubst %.cpp,%,$(wildcard *.cpp))

CXX = g++
CXXFLAGS = -Wall -std=c++11

ifeq ("${DEBUG}","1")
	CXXFLAGS += -g
else
	CXXFLAGS += -DNDEBUG -O2
endif

GEN_VERSION := $(shell bash version.sh version.h.in version.h)

.PHONY: all clean

all: ${TARGET}

clean:
	@rm -fv ${TARGET} ${MODULES:%=%.d} ${MODULES:%=%.o} version.h

${TARGET}: ${MODULES:%=%.o}
	${CXX} ${CXXFLAGS} -o $@ $^

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} -o $@ $<

ifneq ("${MAKECMDGOALS}", "clean")
sinclude ${MODULES:%=%.d}
%.d: %.cpp
	@echo "Parsing dependency for '$<'"
	@${CXX} -MM $< -MT ${@:%.d=%.o} | sed 's,\($*\)\.o[ :]*,\1.o $@: ,g' > $@
endif
