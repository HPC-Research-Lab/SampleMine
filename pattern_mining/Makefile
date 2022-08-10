DEBUG=-O2
OPTIONS=
CXX=g++

HOME_DIR = ..
INCLUDE = -I../ -I$(BOOST_ROOT)
OBJ = join.o match.o filter.o pattern.o patlist.o sampler.o match_asap.o
EXE = sc4_tv.exe sc4_sv.exe sc4_tv_asap.exe sc5_tv.exe sc5_sv.exe sc5_tv_asap.exe sc6_tv.exe sc6_sv.exe sc6_tv_asap.exe sc7_tv.exe sc7_sv.exe sc7_tv_asap.exe fsm4_tv.exe fsm4_sv.exe fsm5_tv.exe fsm5_sv.exe fsm6_tv.exe fsm6_sv.exe 

all: $(EXE)

%.exe: test/%.cpp $(OBJ) bliss/libbliss.a $(wildcard $(HOME_DIR)/db/*.cpp) $(wildcard $(HOME_DIR)/graph/*.cpp) $(HOME_DIR)/util.cpp
	$(CXX)  $(DEBUG) $(OPTIONS) -fopenmp $(INCLUDE) -o $@ $^  -lstdc++fs -ldl -std=c++17 -lstdc++  

$(OBJ): %.o: %.cpp
	$(CXX) $(DEBUG) $(OPTIONS) -fopenmp $(INCLUDE) -c $< -std=c++17 

bliss/libbliss.a:
	$(MAKE) -C bliss 

clean:
	rm *.exe $(EXE) $(OBJ); $(MAKE) -C bliss clean

