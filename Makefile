all:
	g++ huxley.cpp -o huxley -I /usr/local/include/lapackpp/ -llapackpp -DLA_BOUND_CHECK -g  -Wall -pedantic -O6