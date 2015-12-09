TARGET = crabber
MODULES = $(patsubst %.cpp,%,$(wildcard *.cpp))

CXX = g++
CXXFLAGS = -Wall -O2 -std=c++11

.PHONY: all clean

all: ${TARGET}

clean:
	@rm -fv ${TARGET} ${MODULES:%=%.d} ${MODULES:%=%.o}

${TARGET}: ${MODULES:%=%.o}
	${CXX} ${CXXFLAGS} -o $@ $^

%.o: %.cpp
	${CXX} ${CXXFLAGS} -c -o $@ $^

ifneq ("${MAKECMDGOALS}", "clean")
sinclude ${MODULES:%=%.d}
%.d: %.cpp
	@echo "Parsing dependency for '$<'"
	@${CXX} -MM $< -MT ${@:%.d=%.o} | sed 's,\($*\)\.o[ :]*,\1.o $@: ,g' > $@
endif
