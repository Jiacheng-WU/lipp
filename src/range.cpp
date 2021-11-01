#include <lipp.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_set>
#include <random>
#include <unistd.h>
#include <cstring>

using namespace std;

void Example()
{
    LIPP<int64_t, double> lipp;
    lipp.insert(-1, 2);
    lipp.insert(1, 2);
    lipp.insert(2, 3);
    lipp.insert(10, 2);

    int64_t *results = new int64_t[1024];
    {
        int size = lipp.range_query_len(results, 1, 10);
        cout << "Results for range query with given length :" << endl;
        for (int i = 0; i < size; i++)
        {
            cout << results[i] << endl;
        }
    }
}

int main()
{
    Example();

    return 0;
}
