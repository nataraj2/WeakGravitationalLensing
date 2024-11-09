# Compiler
CXX = mpicxx

# Compiler flags
CXXFLAGS = -std=c++14 -Wall -O2 -I. -MMD -MP

# Include directories (if needed)
INCLUDES = -I.

# Source files
SRCS = IO.cpp WL_Functions.cpp main.cpp

# Object directory
OBJDIR = obj

# Object files (stored in obj directory)
OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRCS))

# Dependency files (stored in obj directory)
DEPS = $(OBJS:.o=.d)

# Executable name
EXEC = WL.exe

# Default target
all: $(EXEC)

# Create the object directory
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Linking the executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Compile each .cpp file into .o and store in obj directory
$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

-include $(DEPS)

# Clean up build files
clean:
	rm -rf $(OBJDIR) $(EXEC) *.vtk
