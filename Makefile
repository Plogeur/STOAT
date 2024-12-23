CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -I/usr/include -I/usr/local/include  # Include standard and local paths

# Source files
SRCS = src/matrix.cpp src/arg_parser.cpp src/snarl_parser.cpp src/binary_analysis.cpp src/quantitative_analysis.cpp src/main.cpp

# Output executable name
TARGET = stoat_cxx

# Default rule
all: $(TARGET)

# Rule to build the target
$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

# Clean rule to remove generated files
clean:
	rm -f $(TARGET)
