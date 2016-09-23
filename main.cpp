#include <QCoreApplication>
/*!
  Test barnes hut for a space charge distribution
 */

#include "vectortemplate.h"
#include <iostream>
#include "chargecloud.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    math::Vector<double, 3> vec;

    std::cout << vec << std::endl;

    ChargeCloud<double> cloud;

    cloud.addParticle(1,{1,0,0});
    cloud.addParticle(2,{-1,0,0});

    std::cout << cloud.maxDispersion() << std::endl;

    return a.exec();
}
