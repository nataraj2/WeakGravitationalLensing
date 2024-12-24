# Compiler and flags
CXX = mpicxx
#CXXFLAGS = -std=c++14 -I/usr/local/include -Wall -O2 -I. -MMD -MP
CXXFLAGS = -std=c++14 -I/opt/cray/pe/fftw/3.3.10.3/x86_milan/include -L/opt/cray/pe/fftw/3.3.10.3/x86_milan/lib -Wall -O2 -I. -MMD -MP
LDFLAGS = -L/usr/local/lib -lfftw3

# Source files and target
SOURCES = main.cpp DarkMatter.cpp IO.cpp
HEADERS = DarkMatter.H
TARGET = DM_Init.exe

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Object file rules
main.o: main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

IO.o: IO.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean

