#include "window.h"
#include "approx1.h"
#include "approx2.h"

#include <QKeyEvent>
#include <QPainter>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static const int FUNC_COUNT = 7;

static double func_0(double /*x*/) { return 1.0; }
static double func_1(double x) { return x; }
static double func_2(double x) { return x * x; }
static double func_3(double x) { return x * x * x; }
static double func_4(double x) { return x * x * x * x; }
static double func_5(double x) { return exp(x); }
static double func_6(double x) { return 1.0 / (25.0 * x * x + 1.0); }

static double (*const FUNCS[FUNC_COUNT])(double) = {
    func_0, func_1, func_2, func_3, func_4, func_5, func_6};

static const char *const FUNC_NAMES[FUNC_COUNT] = {
    "f(x)=1",    "f(x)=x",    "f(x)=x^2", "f(x)=x^3",
    "f(x)=x^4", "f(x)=e^x", "f(x)=1/(25x^2+1)"};

static double eval_exact(double t, const Window *w) { return w->f_exact(t); }
static double eval_a1(double t, const Window *w) { return w->eval1(t); }
static double eval_a2(double t, const Window *w) { return w->eval2(t); }
static double eval_e1(double t, const Window *w)
{
    return w->eval1(t) - w->f_exact(t);
}
static double eval_e2(double t, const Window *w)
{
    return w->eval2(t) - w->f_exact(t);
}


static const int DRAW_N = 2000;


Window::Window(QWidget *parent)
    : QWidget(parent), a(-1.0), b(1.0), n(10), k(0), p(0), mode(0), zoom(0),
      x_nodes(nullptr), f_vals(nullptr), coef1(nullptr), coef2(nullptr),
      work(nullptr), alloc_n(0)
{
    setFocusPolicy(Qt::StrongFocus); 
    recompute();
}

Window::~Window() { free_arrays(); }

void Window::free_arrays()
{
    delete[] x_nodes;
    delete[] f_vals;
    delete[] coef1;
    delete[] coef2;
    delete[] work;
    x_nodes = nullptr;
    f_vals = nullptr;
    coef1 = nullptr;
    coef2 = nullptr;
    work = nullptr;
    alloc_n = 0;
}

void Window::allocate(int new_n)
{
    if (new_n == alloc_n) {
        return;
    }
    free_arrays();
    if (new_n < 2) {
        return;
    }
    x_nodes = new double[new_n];
    f_vals = new double[new_n];
    coef1 = new double[4 * (new_n - 1)];
    coef2 = new double[3 * (new_n - 1)];
    work = new double[2 * new_n];
    alloc_n = new_n;
}

void Window::recompute()
{
    int i;
    double max_f, fv;

    allocate(n);
    if (alloc_n < 2) {
        return;
    }

    for (i = 0; i < n; i++) {
        x_nodes[i] = a + (b - a) * i / (n - 1);
        f_vals[i] = FUNCS[k](x_nodes[i]);
    }

    max_f = 0.0;
    for (i = 0; i < n; i++) {
        fv = fabs(f_vals[i]);
        if (fv > max_f) {
            max_f = fv;
        }
    }

    if (p != 0 && n >= 2) {
        f_vals[n / 2] += p * 0.1 * max_f;
    }

    build_approx1(n, x_nodes, f_vals, coef1, work);
    build_approx2(n, x_nodes, f_vals, coef2, work);

    update();
}

int Window::parse_command_line(int argc, char *argv[])
{
    if (argc < 5) {
        return -1;
    }
    if (sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1 ||
        sscanf(argv[3], "%d", &n) != 1 || sscanf(argv[4], "%d", &k) != 1) {
        return -2;
    }
    if (b - a < 1e-9 || n < 2 || k < 0 || k >= FUNC_COUNT) {
        return -3;
    }
    recompute();
    return 0;
}

QSize Window::minimumSizeHint() const { return QSize(200, 200); }

QSize Window::sizeHint() const { return QSize(800, 600); }

void Window::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_0: 
        k = (k + 1) % FUNC_COUNT;
        p = 0;
        recompute();
        break;
    case Qt::Key_1: 
        mode = (mode + 1) % 4;
        update();
        break;
    case Qt::Key_2:
        zoom++;
        update();
        break;
    case Qt::Key_3: 
        zoom--;
        update();
        break;
    case Qt::Key_4:
        n = n * 2;
        p = 0;
        recompute();
        break;
    case Qt::Key_5: 
        if (n / 2 >= 2) {
            n = n / 2;
            p = 0;
            recompute();
        }
        break;
    case Qt::Key_6:
        p++;
        recompute();
        break;
    case Qt::Key_7: 
        p--;
        recompute();
        break;
    default:
        QWidget::keyPressEvent(event);
        break;
    }
}


double Window::f_exact(double t) const { return FUNCS[k](t); }

double Window::eval1(double t) const
{
    if (alloc_n < 2) {
        return 0.0;
    }
    return eval_approx1(t, a, b, n, x_nodes, coef1);
}

double Window::eval2(double t) const
{
    if (alloc_n < 2) {
        return 0.0;
    }
    return eval_approx2(t, a, b, n, x_nodes, coef2);
}


void Window::get_view(double &va, double &vb) const
{
    int i;
    double factor = 1.0;

    if (zoom >= 0) {
        for (i = 0; i < zoom; i++) {
            factor *= 2.0;
        }
        va = a / factor;
        vb = b / factor;
    } else {
        for (i = 0; i < -zoom; i++) {
            factor *= 2.0;
        }
        va = a * factor;
        vb = b * factor;
    }
}

void Window::draw_curve(QPainter *painter, double va, double vb,
                        double y_min, double y_max,
                        double (*eval_fn)(double, const Window *),
                        const QColor &color) const
{
    int i;
    double x0, y0, x1, y1;
    double range_y = y_max - y_min;
    double W = static_cast<double>(width());
    double H = static_cast<double>(height());

    if (range_y < 1e-30) {
        range_y = 1.0;
    }

    QPen pen(color);
    pen.setWidth(0);
    painter->setPen(pen);

    x0 = va;
    y0 = eval_fn(va, this);
    for (i = 1; i <= DRAW_N; i++) {
        double sx0, sy0, sx1, sy1;

        x1 = va + (vb - va) * i / DRAW_N;
        y1 = eval_fn(x1, this);

        sx0 = (x0 - va) / (vb - va) * W;
        sy0 = H - (y0 - y_min) / range_y * H;
        sx1 = (x1 - va) / (vb - va) * W;
        sy1 = H - (y1 - y_min) / range_y * H;

        painter->drawLine(QPointF(sx0, sy0), QPointF(sx1, sy1));

        x0 = x1;
        y0 = y1;
    }
}

void Window::paintEvent(QPaintEvent * /*event*/)
{
    QPainter painter(this);
    double va, vb;
    double y_min, y_max, y_max_abs;
    double t, v;
    int i;
    char info[256];

    get_view(va, vb);

    if (alloc_n < 2) {
        return;
    }

    if (mode == 3) {
        y_min = eval_e1(va, this);
        y_max = y_min;
    } else {
        y_min = eval_exact(va, this);
        y_max = y_min;
    }

    for (i = 0; i <= DRAW_N; i++) {
        t = va + (vb - va) * i / DRAW_N;

        if (mode != 3) {
            v = eval_exact(t, this);
            if (v < y_min) y_min = v;
            if (v > y_max) y_max = v;
        }
        if (mode == 0 || mode == 2) {
            v = eval_a1(t, this);
            if (v < y_min) y_min = v;
            if (v > y_max) y_max = v;
        }
        if (mode == 1 || mode == 2) {
            v = eval_a2(t, this);
            if (v < y_min) y_min = v;
            if (v > y_max) y_max = v;
        }
        if (mode == 3) {
            v = eval_e1(t, this);
            if (v < y_min) y_min = v;
            if (v > y_max) y_max = v;
            v = eval_e2(t, this);
            if (v < y_min) y_min = v;
            if (v > y_max) y_max = v;
        }
    }

    {
        double margin = 0.01 * (y_max - y_min);

        if (margin < 1e-15) {
            margin = 1e-15;
        }
        y_min -= margin;
        y_max += margin;
    }

    y_max_abs = fabs(y_min) > fabs(y_max) ? fabs(y_min) : fabs(y_max);
    printf("view=[%.4g,%.4g]  max|F|=%.6g  k=%d  n=%d  s=%d  p=%d\n",
           va, vb, y_max_abs, k, n, zoom, p);

    {
        double H = static_cast<double>(height());
        double W = static_cast<double>(width());
        double range_y = y_max - y_min;
        double ax = (0.0 - va) / (vb - va) * W;
        double ay = H - (0.0 - y_min) / range_y * H;

        QPen axpen(QColor("lightGray"));
        axpen.setWidth(0);
        painter.setPen(axpen);
        painter.drawLine(QPointF(0, ay), QPointF(W, ay));
        painter.drawLine(QPointF(ax, 0), QPointF(ax, H));
    }

    switch (mode) {
    case 0:
        draw_curve(&painter, va, vb, y_min, y_max, eval_exact, QColor("black"));
        draw_curve(&painter, va, vb, y_min, y_max, eval_a1, QColor("blue"));
        break;
    case 1:
        draw_curve(&painter, va, vb, y_min, y_max, eval_exact, QColor("black"));
        draw_curve(&painter, va, vb, y_min, y_max, eval_a2, QColor("red"));
        break;
    case 2:
        draw_curve(&painter, va, vb, y_min, y_max, eval_exact, QColor("black"));
        draw_curve(&painter, va, vb, y_min, y_max, eval_a1, QColor("blue"));
        draw_curve(&painter, va, vb, y_min, y_max, eval_a2, QColor("red"));
        break;
    case 3:
        draw_curve(&painter, va, vb, y_min, y_max, eval_e1, QColor("blue"));
        draw_curve(&painter, va, vb, y_min, y_max, eval_e2, QColor("red"));
        break;
    default:
        break;
    }

    snprintf(info, sizeof(info),
             "k=%d %s  n=%d  mode=%d  zoom=%d  p=%d  max|F|=%.4g",
             k, FUNC_NAMES[k], n, mode, zoom, p, y_max_abs);
    painter.setPen(QColor("black"));
    painter.drawText(10, 20, QString(info));
    painter.setPen(QColor("blue"));
    painter.drawText(10, 38, QString("blue=Akima(1)"));
    painter.setPen(QColor("red"));
    painter.drawText(10, 54, QString("red=ParabSpline(2)"));
}
