mfcc : main.o hmath.o ini.o mfcc.o sigProcess.o WAVE.o
	gcc -lm -o mfcc main.o hmath.o ini.o mfcc.o sigProcess.o WAVE.o

main.o:main.c ini.h hmath.h mfcc.h sigProcess.h WAVE.h
	gcc -c main.c -lm

hmath.o:hmath.c hmath.h
	gcc -c hmath.c

ini.o:ini.c ini.o
	gcc -c ini.c

mfcc.o:mfcc.c mfcc.h hmath.h
	gcc -c mfcc.c

sigProcess.o:sigProcess.c sigProcess.h hmath.h
	gcc -c sigProcess.c

WAVE.o:WAVE.c WAVE.h hmath.h
	gcc -c WAVE.c

clean:
	rm -f mfcc mfcc *.o
