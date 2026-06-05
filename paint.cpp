#include "paint.h"
#include <math.h>
#include <stdio.h>

static const int DRAW_SAMPLES = 1000;

double draw_function(QPainter *painter, double a, double b,
                     double (*func)(double), const QColor &color)
{
    int i;
    double min_y, max_y, delta_y, scale_x, scale_y;
    double x, y, x_prev, y_prev;
    double max_abs;

    min_y = func(a);
    max_y = min_y;
    for (i = 0; i <= DRAW_SAMPLES; i++) {
        x = a + (b - a) * i / DRAW_SAMPLES;
        y = func(x);
        if (y < min_y) {
            min_y = y;
        }
        if (y > max_y) {
            max_y = y;
        }
    }

    delta_y = 0.01 * (max_y - min_y);
    if (delta_y < 1e-15) {
        delta_y = 1e-15;
    }
    min_y -= delta_y;
    max_y += delta_y;

    max_abs = fabs(min_y) > fabs(max_y) ? fabs(min_y) : fabs(max_y);
    printf("max|f| = %.6g\n", max_abs);

    scale_x = painter->device()->width() / (b - a);
    scale_y = painter->device()->height() / (max_y - min_y);

    QPen pen(color);
    pen.setWidth(0);
    painter->setPen(pen);

    painter->save();
    painter->translate(0.0, painter->device()->height());
    painter->scale(scale_x, -scale_y);
    painter->translate(-a, -min_y);

    x_prev = a;
    y_prev = func(a);
    for (i = 1; i <= DRAW_SAMPLES; i++) {
        x = a + (b - a) * i / DRAW_SAMPLES;
        y = func(x);
        painter->drawLine(QPointF(x_prev, y_prev), QPointF(x, y));
        x_prev = x;
        y_prev = y;
    }

    painter->restore();

    return max_abs;
}
