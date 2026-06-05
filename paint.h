#ifndef PAINT_H
#define PAINT_H

#include <QPainter>

double draw_function(QPainter *painter, double a, double b,
                     double (*func)(double), const QColor &color);

#endif 
