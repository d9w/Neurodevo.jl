CC = mpicxx
SOURCES = core/Axon.cpp core/Soma.cpp core/Environment.cpp
CFLAGS = -g -std=c++14 -Wall -Wextra -Wshadow -pedantic
LFLAGS =
DRAWLFLAGS = -lGL -lglut
OPTFLAGS = -O3
OBJECTS=$(SOURCES:.cpp=.o)

main:$(OBJECTS) run.o
	$(CC) $(CFLAGS) $(OBJECTS) run.o -o run $(LFLAGS) $(OPTFLAGS)

evo:$(OBJECTS) evo.o
	$(CC) $(CFLAGS) $(OBJECTS) evo.o -o evo $(LFLAGS) $(OPTFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@ $(OPTFLAGS)

clean:
	rm -f evo run *.o core/*.o
