P1 = ../BasicAudioToolBox

mfcc : main.o hmath.o ini.o mfcc.o sigProcess.o WAVE.o cnpy.o fileIO.o
	g++ -o mfcc main.o hmath.o ini.o mfcc.o sigProcess.o WAVE.o cnpy.o fileIO.o -fopenmp

main.o:main.c ../inih/ini.h $(P1)/hmath.hpp mfcc.h $(P1)/sigProcess.hpp $(P1)/WAVE.hpp cnpy.hpp
	g++ -c main.c -I ../inih -I $(P1) -fopenmp

mfcc.o:mfcc.c mfcc.h $(P1)/hmath.hpp $(P1)/sigProcess.hpp $(P1)/WAVE.hpp $(P1)/fileIO.hpp cnpy.hpp
	g++ -c mfcc.c -I ../inih -I $(P1)

cnpy.o:cnpy.cpp cnpy.hpp
	g++ -c cnpy.cpp

hmath.o:$(P1)/hmath.cpp $(P1)/hmath.hpp
	g++ -c $(P1)/hmath.cpp

sigProcess.o:$(P1)/sigProcess.cpp $(P1)/sigProcess.hpp $(P1)/hmath.hpp
	g++ -c $(P1)/sigProcess.cpp

WAVE.o:$(P1)/WAVE.cpp $(P1)/WAVE.hpp $(P1)/hmath.hpp
	g++ -c $(P1)/WAVE.cpp

fileIO.o:$(P1)/fileIO.cpp $(P1)/fileIO.hpp
	g++ -c $(P1)/fileIO.cpp

ini.o:../inih/ini.c ../inih/ini.h
	g++ -c ../inih/ini.c

clean:
	rm -f mfcc *.o
