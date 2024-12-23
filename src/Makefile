CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2  # Use C++17, enable warnings, optimize

# Source files
SRCS = matrix.cpp arg_parser.cpp snarl_parser.cpp main.cpp
# Output executable name
TARGET = slink

# Default rule
all: $(TARGET)

# Rule to build the target
$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

# Clean rule to remove generated files
clean:
	rm -f $(TARGET)
