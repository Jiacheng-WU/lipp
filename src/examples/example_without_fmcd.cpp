#include <lipp.h>
#include <iostream>

using namespace std;

int main()
{
    // disable FMCD algorithm
    LIPP<int, int, false> lipp;

    // insert key-values
    lipp.insert(1, 1);
    lipp.insert(2, 3);
    lipp.insert(10, 2);

    cout << "exists(1) = " << (lipp.exists(1) ? "true" : "false") << endl;
    cout << "exists(4) = " << (lipp.exists(4) ? "true" : "false") << endl;

    cout << lipp.at(2) << endl;

    // show tree structure
    lipp.show();

    return 0;
}
