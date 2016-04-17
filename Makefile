CXX = g++
CXXFLAGS = -Wall -g
LDFLAGS =
MAIN = testHeston

SOURCES = testHeston.cpp Option.cpp Vanille.cpp mtrand.cpp MC.cpp MonteCarlo.cpp Heston.cpp BS.cpp
OBJETS = $(SOURCES:.cpp=.o)

all: $(MAIN)


testHeston : $(OBJETS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
.cxx.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean :
	rm -f $(MAIN) $(OBJETS)


