DEBUG=-O2
OPTIONS=
CXX=g++

HOME_DIR = ..
INCLUDE = -I../ -I$(BOOST_ROOT)
OBJ = join.o match.o filter.o pattern.o patlist.o
EXE =  fsm5_nested_sampling.exe fsm5_sampling.exe #mc4_without_sampling.exe mc5_without_sampling.exe mc4_without_sampling_two_vertex.exe mc5_without_sampling_two_vertex.exe fsm4_single_vertex.exe

all: $(EXE)

%.exe: test/%.cpp $(OBJ) bliss/libbliss.a $(wildcard $(HOME_DIR)/db/*.cpp) $(wildcard $(HOME_DIR)/graph/*.cpp) $(HOME_DIR)/util.cpp
	$(CXX)  $(DEBUG) $(OPTIONS) -fopenmp $(INCLUDE) -o $@ $^  -lstdc++fs -ldl -std=c++17 -lstdc++  

$(OBJ): %.o: %.cpp
	$(CXX) $(DEBUG) $(OPTIONS) -fopenmp $(INCLUDE) -c $< -std=c++17 

bliss/libbliss.a:
	$(MAKE) -C bliss 


clean:
	rm *.exe $(EXE) $(OBJ); $(MAKE) -C bliss clean

