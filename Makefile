all:
	g++  tracking.cpp `pkg-config --cflags --libs opencv` && ./a.out