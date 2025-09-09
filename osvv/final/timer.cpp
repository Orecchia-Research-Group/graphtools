#include <chrono>

float timer()
{
    using namespace std::chrono;
    static steady_clock::time_point start = steady_clock::now();

    // elapsed time since first call
    duration<float> elapsed = steady_clock::now() - start;
    return elapsed.count();  // in seconds
}
