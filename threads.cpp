#include <iostream>
#include <vector>
#include <thread>
#include <chrono>

const int SIZE = 100000000; // Size of the array
const int THREADS = std::thread::hardware_concurrency(); // Number of threads

// Function to compute the sum of a segment
void partial_sum(const std::vector<int>& arr, long long& result, int start, int end) {
    result = 0;
    for (int i = start; i < end; ++i) {
        result += arr[i];
    }
}

int main() {
    std::cout << "Available threads: " << THREADS << "\n";
    return 0;
}