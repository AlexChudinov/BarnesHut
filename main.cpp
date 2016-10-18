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
    using matrix1= math::matrix_c<double, 3, 4>;
    using vector_c = math::vector_c<double, 3>;

    matrix m{{1.,2.,3.},{4.,5,6.},{7,8.,9}};
    matrix1 m1{{1.,2.,3.,4.},{5.,6.,7.,8.},{9.,10.,11.,12.}};
    vector_c v{4., 1., 3.};

    //std::cout << math::solve(m,v) << std::endl;
    std::cout << ((v+v)*10.)/2. << std::endl;
    std::cout << (v * (2.*v)) << std::endl;
    std::cout << (m*m)*m << std::endl;
    std::cout << math::transpose(10.*m) << std::endl;
    std::cout << (math::transpose(m1)*m) << std::endl;

    std::vector<vector_c> vv{{1., 2., 3.}, {2., 13., 4.},
                           {4., 15., 6.}, {6., 7., 8.}};

    vector_c v0(0.0);

    std::vector<double> w{1,1,1,1};

    math::mean<vector_c, double, std::vector>(vv,w,v0);

    std::cout << v0 << std::endl;

    std::cout << math::pc1(vv,{1.,1.,1.,1.}) << std::endl;

    /*std::cout << math::det(m) << std::endl;*/

    return a.exec();
}
