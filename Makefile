mfcc : main.o hmath.o ini.o mfcc.o sigProcess.o WAVE.o cnpy.o
	g++ -o mfcc main.o hmath.o ini.o mfcc.o sigProcess.o WAVE.o cnpy.o

main.o:main.c ini.h hmath.h mfcc.h sigProcess.h WAVE.h cnpy.hpp
	g++ -c main.c

hmath.o:hmath.c hmath.h
	g++ -c hmath.c

ini.o:ini.c ini.o
	g++ -c ini.c

mfcc.o:mfcc.c mfcc.h hmath.h
	g++ -c mfcc.c

sigProcess.o:sigProcess.c sigProcess.h hmath.h
	g++ -c sigProcess.c

WAVE.o:WAVE.c WAVE.h hmath.h
	g++ -c WAVE.c

cnpy.o:cnpy.cpp cnpy.hpp
	g++ -c cnpy.cpp

clean:
	rm -f mfcc mfcc *.o
