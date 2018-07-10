CXXFLAGS += -I. $(shell root-config --cflags) -g
CXXFLAGS += -I. Eigen/Dense.h
LDFLAGS += $(shell root-config --libs) -lPhysics -lMatrix -g

PROGRAMS = DistortionClass CalibSCE CalcEField

all:		clean $(PROGRAMS)

$(PROGRAMS):
	@echo '<<compiling' $@'>>'
	@$(CXX) $@.cpp -o $@ $(CXXFLAGS) $(LDFLAGS)
	@rm -rf *.dSYM
clean:	
	rm -f $(PROGRAMS)
