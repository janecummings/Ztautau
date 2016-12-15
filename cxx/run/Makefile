ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBSIM=$(shell root-config --glibs) -L$(ROOTSYS)/lib -lCore -lCint -lHist -lMatrix -ldl -lGpad -lGraf -lTree -lPhysics -lFoam -lTreePlayer 
ROOTLIBS=$(shell root-config --glibs) -L$(ROOTSYS)/lib -lCore -lCint -lHist -lMatrix -ldl -lGpad -lGraf -lTree -lPhysics -lFoam -lTreePlayer 
ROOTCINT=$(ROOTSYS)/bin/rootcint
CXX=g++ -fPIC 
CXXFLAGS=-I$(ROOTSYS)/include 
LDLFLAGS=$(ROOTLIBS)

TOP := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
TAUSPINNER=/afs/cern.ch/user/c/cummings/TauSpinner/TauSpinnerTool/tags/TauSpinnerTool-00-00-04
EXTERN=extern
HEPMC=$(TAUSPINNER)/$(EXTERN)/src/HEPMC
TAUOLA=$(TAUSPINNER)/$(EXTERN)/src/TAUOLA
LHAPDF=$(TAUSPINNER)/pdfs
ROOTLIBS += -L$(TAUSPINNER)/lib -lTauSpinnerTool -L$(TAUSPINNER)/$(EXTERN)/lib -lTauolaFortran -lTauolaCxxInterface -lHepMC -lLHAPDF -lTauolaTauSpinner
CXXFLAGS += -I$(TAUSPINNER)/$(EXTERN)/include -I$(TAUSPINNER)/$(EXTERN)/include/Tauola -I$(TAUSPINNER)/src

HEADERS=SpinTree.h
MODULES=Dict.o SpinTree.o tauspin.o
SIMPLE=Dict.o SpinTree.o spin.o
PYTHIA=Dict.o PythiaTree.o recopythia.o
COEF=Dict.o SpinTree.o makeCoef.o
SPINNER=Dict.o SpinTree.o spinner.o
REWEIGHT=Dict.o SpinTree.o reweight.o
COSTHETA=Dict.o SpinTree.o costheta.o
ZPT=Dict.o SpinTree.o Zptspin.o


Dict.cxx: $(HEADERS) LinkDef.h
	@echo "Generating Dictionary ...."
	rootcint -f $@ -c $^

simple: $(SIMPLE)
	$(CXX) $(CXXFLAGS) $(ROOTLIBSIM) $(SIMPLE) -o spin

tauspin: $(MODULES)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) -Wl,-rpath,$(ROOTSYS)/lib -Wl,-rpath,$(TAUSPINNER)/lib -Wl,-rpath,$(TAUSPINNER)/$(EXTERN)/lib $(MODULES) -o tauspin

Zptspin: $(ZPT)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) -Wl,-rpath,$(ROOTSYS)/lib -Wl,-rpath,$(TAUSPINNER)/lib -Wl,-rpath,$(TAUSPINNER)/$(EXTERN)/lib $(ZPT) -o Zptspin

pythiaspin: $(PYTHIA)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) -Wl,-rpath,$(ROOTSYS)/lib -Wl,-rpath,$(TAUSPINNER)/lib -Wl,-rpath,$(TAUSPINNER)/$(EXTERN)/lib $(PYTHIA) -o pythiaspin

recopythia: $(PYTHIA)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) -Wl,-rpath,$(ROOTSYS)/lib -Wl,-rpath,$(TAUSPINNER)/lib -Wl,-rpath,$(TAUSPINNER)/$(EXTERN)/lib $(PYTHIA) -o recopythia

coef: $(COEF)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) -Wl,-rpath,$(ROOTSYS)/lib -Wl,-rpath,$(TAUSPINNER)/lib -Wl,-rpath,$(TAUSPINNER)/$(EXTERN)/lib $(COEF) -o makecoef

spinner: $(SPINNER)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) -Wl,-rpath,$(ROOTSYS)/lib -Wl,-rpath,$(TAUSPINNER)/lib -Wl,-rpath,$(TAUSPINNER)/$(EXTERN)/lib $(SPINNER) -o spinner

reweight: $(REWEIGHT)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) -Wl,-rpath,$(ROOTSYS)/lib -Wl,-rpath,$(TAUSPINNER)/lib -Wl,-rpath,$(TAUSPINNER)/$(EXTERN)/lib $(REWEIGHT) -o reweight

costheta: $(COSTHETA)
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) -Wl,-rpath,$(ROOTSYS)/lib -Wl,-rpath,$(TAUSPINNER)/lib -Wl,-rpath,$(TAUSPINNER)/$(EXTERN)/lib $(COSTHETA) -o costheta

%.o: %.cxx %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -rf *.o 
