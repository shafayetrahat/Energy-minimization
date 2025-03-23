CXX = g++
CXXFLAGS = -g -Wall -O2
LDFLAGS =
SRCDIR = src
BINDIR = bin


all: lj_2d lj_3d linked_cell_2d

lj_2d: $(SRCDIR)/lj_2d.cpp
	mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $< $(LDFLAGS)

lj_3d: $(SRCDIR)/lj_3d.cpp
	mkdir -p $(BINDIR) 
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $< $(LDFLAGS)

linked_cell_2d: $(SRCDIR)/linked_cell_2d.cpp
	mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $< $(LDFLAGS)
clean:
	rm -rf $(BINDIR)/*

cleanall:
	rm -rf run/* bin/*

test:
	mkdir -p run
	./bin/lj_2d
	./bin/lj_3d
