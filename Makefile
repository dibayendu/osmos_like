#Makefile: Sun Jul 22 00:24:15 EST 2012 pknowles */

DEBUG = -g
#OPTIMISE = -O3

CFLAGS = `pkg-config gdk-pixbuf-2.0 --cflags` `sdl-config --cflags` $(DEBUG) $(OPTIMISE) -std=c99
LDFLAGS = `pkg-config gdk-pixbuf-2.0 --libs` `sdl-config --libs` -lGL -lGLU -lglut -lm -lz
EXE = assignment_3_application
ZIP = 3224807_assignment3.zip
OBJECTS = main.o sdlbase.o arena.o particle.o collision.o utility.o texture.o

#default target
all: $(EXE)

#executable
$(EXE) : $(OBJECTS)
	gcc -o $@ $(LDFLAGS) $(OBJECTS)

#general object file rule
%.o : %.c
	gcc -c -o $@ $(CFLAGS) $<

#additional dependencies
sdlbase.o : sdlbase.h
utility.o: utility.h
main.o : sdlbase.h
arena.o: utility.h
particle.o: utility.h
collision.o: utility.h particle.h
texture.o: texture.h

#clean (-f stops error if files don't exist)
clean:
	rm -f $(EXE) \
	      $(OBJECTS) \
		  $(ZIP)

archive:
	zip 3224807_assignment3.zip *.h *.c *.txt Makefile -r images/

square:
	./assignment_3_application 1 # 1 is for starting a square arena

circle:
	./assignment_3_application 0 # 0 is for starting a circle arena
