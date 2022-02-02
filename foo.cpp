#include <iostream>
using namespace std;

extern "C" void cfoo(const unsigned long *indatav)
{
    for (int i = 0; i < 64; ++i){
        cout << "c++ vals --> " << indatav[i] << '\n';
    }
}