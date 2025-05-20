SRCDIR    = src
INCDIR    = include
OBJDIR    = build
CXX       = g++
CPPFLAGS  = -I./${INCDIR}
CXXFLAGS  = -O2 -Wall
LDFLAGS   = -lyaml-cpp -lcln -lginac


OBJS      = ${OBJDIR}/amflow.cpp.o    \
			${OBJDIR}/apart.cpp.o     \
			${OBJDIR}/family.cpp.o    \
			${OBJDIR}/integral.cpp.o  \
			${OBJDIR}/ibp.cpp.o       \
			${OBJDIR}/kira.cpp.o      \
			${OBJDIR}/reduction.cpp.o \
			${OBJDIR}/boundary.cpp.o  \
			${OBJDIR}/main.cpp.o


all: pre amflow

amflow: ${OBJS}
	${CXX} ${OBJS} -o amflow ${LDFLAGS}

${OBJS}: ${OBJDIR}/%.cpp.o: ${SRCDIR}/%.cpp
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $^ -o $@

.PHONY: pre
pre:
	mkdir -p ${OBJDIR}

.PHONY: clean
clean:
	rm -rf ${OBJDIR}
	rm -f amflow


