QT += core
QT -= gui

CONFIG += c++11

TARGET = BarnesHut
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp

HEADERS += \
    chargecloud.h \
    vectortemplate.h \
    mathutility.h \
    matrixtemplate.h
