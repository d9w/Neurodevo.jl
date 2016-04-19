CC = mpicxx
SOURCES = core/Axon.cpp core/Soma.cpp core/DNA.cpp core/ANN.cpp problems/Robot.cpp problems/Forage.cpp
CFLAGS = -g -std=c++14 -Wall -Wextra -Wshadow -pedantic
LFLAGS =
DRAWLFLAGS = -lGL -lGLEW -lGLU -lglfw
OPTFLAGS = -O3
OBJECTS=$(SOURCES:.cpp=.o)

main:$(OBJECTS) run.o
	$(CC) $(CFLAGS) $(OBJECTS) run.o -o run $(LFLAGS) $(OPTFLAGS)

evo:$(OBJECTS) evo.o
	$(CC) $(CFLAGS) $(OBJECTS) evo.o -o evo $(LFLAGS) $(OPTFLAGS)

draw:$(OBJECTS) core/Viewer.o draw.o
	$(CC) $(CFLAGS) $(OBJECTS) core/Viewer.o draw.o -o draw $(LFLAGS) $(DRAWLFLAGS) $(OPTFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@ $(OPTFLAGS)

clean:
	rm -f evo run draw *.o core/*.o problems/*.o
