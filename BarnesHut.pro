QT += core
QT -= gui

CONFIG += c++11
#QMAKE_CXXFLAGS -= -O2
#QMAKE_CXXFLAGS += -O3

TARGET = BarnesHut
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp

HEADERS += \
    chargecloud.h \
    vectortemplate.h \
    mathutility.h \
    matrixtemplate.h \
    binarytree.h \
    array_operations.h

INCLUDEPATH += D:/boost_1_62_0
