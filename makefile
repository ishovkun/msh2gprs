CXX = g++
CXXFLAGS = -pipe -O2 -m64 -std=c++0x

PREDCXXFLAGS = -O2

SWITCHES = -DSELF_CHECK

customtetgen:	main.cpp simdata.cpp transes.cpp femout.cpp element.cpp renum.cpp Polyhedra.cpp Collision.cpp Vector.cpp
	$(CXX) $(CXXFLAGS) $(SWITCHES) -o ./bin/customgmsh main.cpp simdata.cpp transes.cpp femout.cpp element.cpp renum.cpp Polyhedra.cpp Collision.cpp Vector.cpp

clean:
	rm -r *.o customgmsh *~
