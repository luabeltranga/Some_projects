GCC=g++
FLAGS= -fsanitize=leak -fsanitize=address -fsanitize=undefined -Wall

wave1d:wave1d.cpp
	$(GCC) $(FLAGS) wave1d.cpp
wave2d:wave2d.cpp
	$(GCC) $(FLAGS) wave2d.cpp

clean:
	rm -rf *.out *.gif *.txt *.dat
