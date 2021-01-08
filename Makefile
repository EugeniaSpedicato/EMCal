LDFLAGS=$(shell root-config --libs)
INCLUDE=$(shell root-config --incdir)
CXXFLAGS=$(shell root-config --cflags)
CXX=$(shell root-config --cxx)
main: main.cc ECAL.cc EMShower.cc RadialInterval.cc GammaFunctionGenerator.cc BaseNumericalRandomGenerator.cc IncGamma.cc
	${CXX} -std=c++1y ${CXXFLGS} -I${INCLUDE} ${LDFLAGS} -o main main.cc ECAL.cc EMShower.cc RadialInterval.cc GammaFunctionGenerator.cc BaseNumericalRandomGenerator.cc IncGamma.cc
