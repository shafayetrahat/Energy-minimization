CXX = g++
CXXFLAGS = -g -std=c++17 -Wall -O3
LDFLAGS =
SRCDIR = src
BINDIR = bin
RUNDIR = run

all: lj_2d lj_3d  lj_2d_variable_box lj_3d_variable_box linked_cell_2d linked_cell_2d_variable_box linked_cell_3d linked_cell_3d_variable_box

lj_2d: $(SRCDIR)/lj_2d.cpp
	mkdir -p $(BINDIR) $(RUNDIR) 
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $< $(LDFLAGS)

lj_2d_variable_box: $(SRCDIR)/lj_2d_variable_box.cpp
	mkdir -p $(BINDIR) $(RUNDIR)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $< $(LDFLAGS)

lj_3d: $(SRCDIR)/lj_3d.cpp
	mkdir -p $(BINDIR) $(RUNDIR)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $< $(LDFLAGS)

lj_3d_variable_box: $(SRCDIR)/lj_3d_variable_box.cpp
	mkdir -p $(BINDIR) $(RUNDIR)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $< $(LDFLAGS)

linked_cell_2d: $(SRCDIR)/linked_cell_2d.cpp
	mkdir -p $(BINDIR) $(RUNDIR)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $< $(LDFLAGS)

linked_cell_2d_variable_box: $(SRCDIR)/linked_cell_2d_variable_box.cpp
	mkdir -p $(BINDIR) $(RUNDIR)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $< $(LDFLAGS)

linked_cell_3d: $(SRCDIR)/linked_cell_3d.cpp
	mkdir -p $(BINDIR) $(RUNDIR)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $< $(LDFLAGS)

linked_cell_3d_variable_box: $(SRCDIR)/linked_cell_3d_variable_box.cpp
	mkdir -p $(BINDIR) $(RUNDIR)
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $< $(LDFLAGS)
	
clean:
	rm -rf $(BINDIR)/*

cleanall:
	rm -rf $(RUNDIR)/* $(BINDIR)/*

test:
	./bin/lj_2d
	./bin/lj_3d
	./bin/lj_2d_variable_box
	./bin/lj_3d_variable_box
	./bin/linked_cell_2d
	./bin/linked_cell_2d_variable_box
	./bin/linked_cell_3d
	./bin/linked_cell_3d_variable_box