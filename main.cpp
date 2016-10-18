#include <QCoreApplication>
/*!
  Test barnes hut for a space charge distribution
 */

#include "chargecloud.h"
#include <iostream>
//#include <ctime>
#include "vectortemplate.h"
#include "matrixtemplate.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    ChargeCloud<double> cloud({1.,1.},{{1.,2.,3.},{-1.,3.,4.}});

    using barnes_hut_tree = ChargeCloud<double>::barnes_hut_tree;

    barnes_hut_tree tree = cloud.create_barnes_hut_tree(0.1);

    return a.exec();
}
