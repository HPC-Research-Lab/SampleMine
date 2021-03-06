DEBUG=-O2
OPTIONS=
CXX=g++

HOME_DIR = ..
INCLUDE = -I../ -I$(BOOST_ROOT)
OBJ = join.o match.o filter.o pattern.o patlist.o
EXE = fsm5_s_sv.exe fsm5_s_ms.exe #fsm5_s.exe fsm5_s_ms.exe #fsm9_s_sv.exe fsm9_s.exe #fsm5_s.exe #fsm7_s.exe fsm7_s_sv.exe fsm4_s.exe fsm4_s_sv.exe #fsm7_s.exe #fsm4_s.exe fsm5_s.exe #fsm7_ns_sv.exe fsm4_ns_sv.exe fsm5_ns_sv.exe fsm6_ns_sv.exe #fsm_s.exe fsm_ns.exe automine_fsm.exe automine.exe #mc.exe fsm.exe #fsm5.exe fsm6.exe fsm.exe #automine.exe  #motif4.exe motif5.exe match_test.exe #motif5.exe #motif7.exe motif7_direct.exe motif5_direct.exe #motif5_direct.exe motif4_direct.exe motif3_direct.exe motif6.exe motif6_direct.exe spanning_tree_test.exe all_pattern_test.exe #motif5_direct.exe spanning_tree_test.exe motif9.exe 
#INCLUDE += -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/include

all: $(EXE)

%.exe: test/%.cpp $(OBJ) bliss/libbliss.a $(wildcard $(HOME_DIR)/db/*.cpp) $(wildcard $(HOME_DIR)/graph/*.cpp) $(HOME_DIR)/util.cpp
	$(CXX)  $(DEBUG) $(OPTIONS) -fopenmp $(INCLUDE) -o $@ $^  -lstdc++fs -ldl -std=c++17 -lstdc++  

$(OBJ): %.o: %.cpp
	$(CXX) $(DEBUG) $(OPTIONS) -fopenmp $(INCLUDE) -c $< -std=c++17 

bliss/libbliss.a:
	$(MAKE) -C bliss 


clean:
	rm *.exe $(EXE) $(OBJ); $(MAKE) -C bliss clean

