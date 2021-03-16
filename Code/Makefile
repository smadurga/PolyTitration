# Compilation of the MC_Polymer program

main.x : main.o lectura_input.o angle_rotations.o inicialitzacio.o escritura_output.o energy.o analysis.o Monte_Carlo.o random.o
	mpicxx  -o main.x main.o lectura_input.o angle_rotations.o inicialitzacio.o escritura_output.o energy.o analysis.o Monte_Carlo.o random.o -lm

main.o : main.cpp global.h lib.h
	mpicxx  -c main.cpp -lm

lectura_input.o : lectura_input.cpp global.h lib.h
	mpicxx  -c lectura_input.cpp -lm

angle_rotations.o : angle_rotations.cpp global.h lib.h
	mpicxx  -c angle_rotations.cpp -lm

inicialitzacio.o : inicialitzacio.cpp global.h lib.h
	mpicxx  -c inicialitzacio.cpp -lm

escritura_output.o: escritura_output.cpp global.h lib.h
	mpicxx  -c escritura_output.cpp -lm

energy.o: energy.cpp global.h lib.h 
	mpicxx   -c energy.cpp -lm

analysis.o: analysis.cpp global.h lib.h
	mpicxx  -c analysis.cpp -lm

Monte_Carlo.o: Monte_Carlo.cpp global.h lib.h 
	mpicxx  -c Monte_Carlo.cpp -lm

random.o: random.cpp global.h lib.h
	mpicxx  -c random.cpp -lm

.PHONY : run
run: main.x 
	mpirun -np 1 main.x

.PHONY : clean
clean: 
	rm *.o
