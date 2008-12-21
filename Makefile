all:
	g++ particles.cpp -o particles -I /usr/local/include/lapackpp/ -llapackpp -DLA_BOUND_CHECK -g  -Wall -pedantic -O6