DEBUG=-O2
OPTIONS=
CXX=g++

HOME_DIR = ..
INCLUDE = -I../ -I$(BOOST_ROOT)
OBJ = join.o match.o filter.o pattern.o patlist.o sampler.o match_asap.o match2.o 
EXE = automine_fsm.exe automine_sc.exe fsm4_sv.exe fsm5_tv.exe q1_size4.exe q2_size7.exe q4_size4.exe q5_size7.exe sc4_tv.exe sc5_tv.exe sc6_tv.exe sc7_tv.exe fsm4_tv.exe fsm6_sv.exe q1_size7.exe q1_size5.exe q2_size5.exe q3_size4.exe q4_size7.exe sc4_sv.exe sc5_sv.exe sc6_sv.exe sc7_sv.exe fsm5_sv.exe fsm6_tv.exe q2_size4.exe q3_size7.exe q5_size4.exe sc4_tv_asap.exe sc5_tv_asap.exe sc6_tv_asap.exe sc7_tv_asap.exe table8_q1_size4_FR.exe table8_q1_size4_UK.exe table8_q1_size5_MI.exe table8_q1_size5_OK.exe table7.exe table8_q2_size4_FR.exe table8_q2_size4_UK.exe table8_q2_size5_OK.exe

all: $(EXE)

%.exe: PACT_Test/%.cpp $(OBJ) bliss/libbliss.a $(wildcard $(HOME_DIR)/db/*.cpp) $(wildcard $(HOME_DIR)/graph/*.cpp) $(HOME_DIR)/util.cpp
	$(CXX)  $(DEBUG) $(OPTIONS) -fopenmp $(INCLUDE) -o $@ $^  -lstdc++fs -ldl -std=c++17 -lstdc++  

$(OBJ): %.o: %.cpp
	$(CXX) $(DEBUG) $(OPTIONS) -fopenmp $(INCLUDE) -c $< -std=c++17 

bliss/libbliss.a:
	$(MAKE) -C bliss 

clean:
	rm *.exe $(EXE) $(OBJ); $(MAKE) -C bliss clean

