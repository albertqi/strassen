strassen: strassen.cpp
	g++ -std=c++17 -O3 -o strassen strassen.cpp

clean:
	rm -f strassen strassen.o