CC = mpicxx
SOURCES = Protein.cpp GRN.cpp Axon.cpp Soma.cpp Environment.cpp MNIST.cpp
CFLAGS = -std=c++11 -Wall -g
LFLAGS = -lpthread -lrt
SNAP =
#SNAP = ../libraries/snap/snap-core/Snap.o
#ISNAP = -I ../libraries/snap/snap-core -I ../libraries/snap/snap-adv -I ../libraries/snap/glib-core/ -I ../libraries/snap/snap-exp
DRAWLFLAGS = -lGL -lglut
OPTFLAGS = -O2
OBJECTS=$(SOURCES:.cpp=.o)

main:$(OBJECTS) Robot.o Spikes.o GA.o main.o
	$(CC) $(CFLAGS) $(OBJECTS) $(SNAP) Robot.o Spikes.o GA.o main.o -o main $(LFLAGS) $(OPTFLAGS)

test:$(OBJECTS) main_test.o
	$(CC) $(CFLAGS) $(OBJECTS) $(SNAP) main_test.o -o test $(LFLAGS) $(OPTFLAGS)

draw:$(OBJECTS) Robot.o Spikes.o main_draw.o
	$(CC) $(CFLAGS) $(OBJECTS) $(SNAP) Robot.o Spikes.o main_draw.o -o draw $(LFLAGS) $(DRAWLFLAGS) $(OPTFLAGS)

robot:$(OBJECTS) Robot.o Spikes.o main_robot.o
	$(CC) $(CFLAGS) $(OBJECTS) $(SNAP) Robot.o Spikes.o main_robot.o -o robot $(LFLAGS) $(DRAWLFLAGS) $(OPTFLAGS)

random:$(OBJECTS) random_experiment.o
	$(CC) $(CFLAGS) $(OBJECTS) $(SNAP) random_experiment.o -o random $(LFLAGS) $(OPTFLAGS)

show_grn: Protein.o GRN.o show_grn.o
	$(CC) $(CFLAGS) Protein.o GRN.o show_grn.o -o show_grn $(LFLAGS) $(OPTFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@ -w $(ISNAP) $(OPTFLAGS)

clean:
	rm -f *.o main test draw
