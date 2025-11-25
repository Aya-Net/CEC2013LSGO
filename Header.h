#ifndef _HEADER_H
#define _HEADER_H


#include <iostream>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <sstream>

// benchmark set header files
#include "benchmarks/F1.h"
#include "benchmarks/F2.h"
#include "benchmarks/F3.h"
#include "benchmarks/F4.h"
#include "benchmarks/F5.h"
#include "benchmarks/F6.h"
#include "benchmarks/F7.h"
#include "benchmarks/F8.h"
#include "benchmarks/F9.h"
#include "benchmarks/F10.h"
#include "benchmarks/F11.h"
#include "benchmarks/F12.h"
#include "benchmarks/F13.h"
#include "benchmarks/F14.h"
#include "benchmarks/F15.h"

using namespace std;

Benchmarks* generateFuncObj(int funcID);

#endif
