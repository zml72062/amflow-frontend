SRCDIR    = src
INCDIR    = include
OBJDIR    = build
CXX       = g++
CPPFLAGS  = -I./${INCDIR} -I./linear-ode/${INCDIR}
CXXFLAGS  = -O2 -Wall
LDFLAGS   = -lyaml-cpp -lcln -lginac -lflint


OBJS      = ${OBJDIR}/amflow.cpp.o    \
			${OBJDIR}/apart.cpp.o     \
			${OBJDIR}/family.cpp.o    \
			${OBJDIR}/integral.cpp.o  \
			${OBJDIR}/ibp.cpp.o       \
			${OBJDIR}/kira.cpp.o      \
			${OBJDIR}/reduction.cpp.o \
			${OBJDIR}/boundary.cpp.o  \
			${OBJDIR}/border.cpp.o    \
			${OBJDIR}/singlesetup.cpp.o\
			${OBJDIR}/solve.cpp.o     \
			${OBJDIR}/main.cpp.o

DEOBJS	  = linear-ode/${OBJDIR}/interface.cpp.o \
			linear-ode/${OBJDIR}/base.cpp.o      \
			linear-ode/${OBJDIR}/jordan.cpp.o    \
			linear-ode/${OBJDIR}/ratsolver.cpp.o \
			linear-ode/${OBJDIR}/symdiffeq.cpp.o \
			linear-ode/${OBJDIR}/symsolver.cpp.o

DESRCS	  = linear-ode/${SRCDIR}/interface.cpp   \
			linear-ode/${SRCDIR}/base.cpp        \
			linear-ode/${SRCDIR}/jordan.cpp      \
			linear-ode/${SRCDIR}/ratsolver.cpp   \
			linear-ode/${SRCDIR}/symdiffeq.cpp   \
			linear-ode/${SRCDIR}/symsolver.cpp  

all: pre amflow

amflow: ${OBJS} ${DESRCS}
	cd linear-ode; make; cd ..
	${CXX} ${DEOBJS} ${OBJS} -o amflow ${LDFLAGS}

${OBJS}: ${OBJDIR}/%.cpp.o: ${SRCDIR}/%.cpp
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $^ -o $@

.PHONY: pre
pre:
	mkdir -p ${OBJDIR}

.PHONY: clean
clean:
	rm -rf ${OBJDIR}
	rm -rf linear-ode/${OBJDIR}
	rm -f amflow


