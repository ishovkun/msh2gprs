CXX = g++
CXXFLAGS = -pipe -O2 -m64 -std=c++0x

PREDCXXFLAGS = -O2

SWITCHES = -DSELF_CHECK

customtetgen:	main.cpp simdata.cpp transes.cpp femout.cpp element.cpp renum.cpp Plane.cpp Polygon.cpp
	$(CXX) $(CXXFLAGS) $(SWITCHES) -o ./bin/customgmsh main.cpp simdata.cpp transes.cpp femout.cpp element.cpp renum.cpp Plane.cpp Polygon.cpp

clean:
	rm -r *.o customgmsh *~
