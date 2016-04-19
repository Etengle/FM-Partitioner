#ifndef __NET_H__
#define __NET_H__
#include <iostream>
#include <vector>
#include <string>
#include "cell.h"
using namespace std;

class Net {

public:
    Net(string str) : A(0), B(0), name(str){
    } 
    ~Net() {}
    string name;
    int A, B;
    vector <int> cellList;
};

#endif
