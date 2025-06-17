CXX      := g++
CXXFLAGS := -std=c++17 -O3 -fopenmp -Wall -I./liblbfgs

SRCS     := $(wildcard *.cpp)
OBJS     := $(SRCS:.cpp=.o)

.PHONY: all lab2 clean

all: lab2

lab2: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -llbfgs

run: lab2
	@./lab2

video:
	ffmpeg -framerate 25 \
	       -start_number 0 \
	       -i frames/animation%d.png \
	       -c:v libx264 -pix_fmt yuv420p \
	       fluid_simulation.mp4



%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f lab2 $(OBJS)
