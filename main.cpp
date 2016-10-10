#include <QCoreApplication>
/*!
  Test barnes hut for a space charge distribution
 */

#include "vectortemplate.h"
#include "matrixtemplate.h"
#include <iostream>
#include <ctime>


int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    math::vector_c<double,3> v1{1,2,3}, v2{4,10,6}, v3{10,10,10}, v0(0.0);
    std::vector<math::vector_c<double,3>> vv{v1, v2, v3};
    math::mean(vv, std::vector<double>{1, 1, 1}, v0);
    vv[0] -= v0; vv[1] -= v0; vv[2] -= v0;

    math::matrix_c<double,3,3> m1(1.0);
    math::matrix_c<double,3,3> m2(1.0);
    m2[0][2] = 3; m2[0][1] = 2;
    std::cout << m2 << std::endl;
    std::cout << math::transpose(m2) << std::endl;
    //std::cout << (v1*m1) << std::endl;
    //std::cout << (v1/=10.) << std::endl;
    std::cout << v1*m1 << std::endl;

    math::matrix_c_op<double, 2, 3>().component_remove(m2, v1);

    math::vector_c<double, 3> vpc = math::pc1(vv,{1,1,1});

    std::cout << vpc << std::endl;

    return a.exec();
}
