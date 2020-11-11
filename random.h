#ifndef RAND_H
#define RAND_H

#include <math.h>
#include <stdlib.h>

//产生double类型随机数
inline double RandDouble()
{
    double r = static_cast<double>(rand());
    return r / RAND_MAX;
}
//反正是产生随机数
inline double RandNormal()
{
    double x1, x2, w;
    do{
        x1 = 2.0 * RandDouble() - 1.0;
        x2 = 2.0 * RandDouble() - 1.0;
        w = x1 * x1 + x2 * x2;
    }while( w >= 1.0 || w == 0.0);

    w = sqrt((-2.0 * log(w))/w);
    return x1 * w;
}

#endif // random.h