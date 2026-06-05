QMAKE_CXXFLAGS += -Werror
CONFIG += c++14
HEADERS       = window.h \
                approx1.h \
                approx2.h \
                paint.h
SOURCES       = main.cpp \
                window.cpp \
                approx1.cpp \
                approx2.cpp \
                paint.cpp
QT += widgets


