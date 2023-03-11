strassen: strassen.cpp
	g++ -std=c++11 -O3 -o strassen strassen.cpp

clean:
	rm -f strassen strassen.o