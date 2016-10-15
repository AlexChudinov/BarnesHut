#include <QCoreApplication>
/*!
  Test barnes hut for a space charge distribution
 */

//#include "chargecloud.h"
#include <iostream>
//#include <ctime>
#include "vectortemplate.h"
#include "matrixtemplate.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    /*ChargeCloud<double> cloud({1.,1.},{{1.,2.,3.},{-1.,3.,4.}});

    using barnes_hut_tree = ChargeCloud<double>::barnes_hut_tree;

    barnes_hut_tree tree = cloud.create_barnes_hut_tree(0.1);*/

    using matrix = math::matrix_c<double, 3, 3>;
    using vector = math::vector_c<double, 3>;

    matrix m{{1.,2.,3.},{4.,5,6.},{7,8.,9}};
    vector v{4., 1., 3.};

    //std::cout << math::solve(m,v) << std::endl;
    std::cout << ((v+v)*10.)/2. << std::endl;
    std::cout << (v * (2.*v)) << std::endl;
    std::cout << (m*v) << std::endl;

    return a.exec();
}
