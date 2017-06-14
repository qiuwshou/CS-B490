all: CImg.h p1.cpp
	g++ p1.cpp -o p1 -lX11 -lpthread

mac: CImg.h p1.cpp
	g++ p1.cpp -o p1 -I/usr/X11R6/include -L/usr/X11R6/lib -lX11 -lpthread

clean:
	rm p1
