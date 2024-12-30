ifndef BOOST_ROOT
$(error BOOST_ROOT is not set)
endif
all: reverse_mode_automatic_differentiation.exe
#CXX=g++
CXXFLAGS=-std=c++17 -DNDEBUG -O3 -march=native -flto=auto -isystem $(BOOST_ROOT)/include -MMD -MP
OBJECTS=reverse_mode_automatic_differentiation.o 

DEPS=$(OBJECTS:.o=.d)
%.o %.d: %.cpp $(BOOST_ROOT)/include
	$(CXX) -c $< -o $*.o $(CXXFLAGS)

reverse_mode_automatic_differentiation.exe: reverse_mode_automatic_differentiation.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	@find . -type f -name "*.d"|xargs rm -f
	@find . -type f -name "*.o"|xargs rm -f
	@find . -type f -name "*.exe"|xargs rm -f

$(BOOST_ROOT)/include:
	@echo \$$\(BOOST_ROOT\)/include does not exist!
	@exit 1

-include $(DEPS)
