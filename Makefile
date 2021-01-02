LDFLAGS=-L/Applications/root_v6.20.00/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -stdlib=libc++ -lpthread -lm -ldl -lRooFit -lRooFitCore
INCLUDE=/Applications/root_v6.20.00/include
CXXFLAGS=-stdlib=libc++ -pthread -std=c++11 -m64 -I/Applications/root_v6.20.00/include
CXX=c++

main: main.cc ECAL.cc EMShower.cc RadialInterval.cc GammaFunctionGenerator.cc BaseNumericalRandomGenerator.cc IncGamma.cc
	${CXX} -std=c++1y ${CXXFLGS} -I${INCLUDE} ${LDFLAGS} -o main main.cc ECAL.cc EMShower.cc RadialInterval.cc GammaFunctionGenerator.cc BaseNumericalRandomGenerator.cc IncGamma.cc