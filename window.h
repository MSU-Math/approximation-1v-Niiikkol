#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>

class Window : public QWidget
{
    Q_OBJECT

  public:
    explicit Window(QWidget *parent = nullptr);
    ~Window();

    int parse_command_line(int argc, char *argv[]);

    QSize minimumSizeHint() const;
    QSize sizeHint() const;

    double f_exact(double t) const;
    double eval1(double t) const;
    double eval2(double t) const;

  protected:
  protected:
    bool event(QEvent *event) override;
    void paintEvent(QPaintEvent *event);
    void keyPressEvent(QKeyEvent *event);

  private:
    double a; /* original left endpoint */
    double b; /* original right endpoint */
    int n;    /* number of interpolation nodes */
    int k;    /* function index 0..6 */
    int p;    /* perturbation level */
    int mode; /* display mode 0..3 */
    int zoom; /* zoom exponent s (view = [a/2^s, b/2^s]) */

    double *x_nodes; /* interpolation nodes [n] */
    double *f_vals;  /* (possibly perturbed) function values [n] */
    double *coef1;   /* Akima coefficients [4*(n-1)] */
    double *coef2;   /* parabolic coefficients [3*(n-1)] */
    double *work;    /* work array [2*n] */

    int alloc_n; /* currently allocated size */

    void recompute();
    void allocate(int new_n);
    void free_arrays();

    void get_view(double &va, double &vb) const;

    void draw_curve(QPainter *painter, double va, double vb,
                    double y_min, double y_max,
                    double (*eval_fn)(double, const Window *),
                    const QColor &color) const;
};

#endif /* WINDOW_H */
