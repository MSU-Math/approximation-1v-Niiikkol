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
    void paintEvent(QPaintEvent *event);
    void keyPressEvent(QKeyEvent *event);

  private:
    double a; 
    double b; 
    int n;    
    int k;    
    int p;    
    int mode; 
    int zoom; 

    double *x_nodes; 
    double *f_vals;  
    double *coef1;   
    double *coef2;   
    double *work;    

    int alloc_n; 

    void recompute();
    void allocate(int new_n);
    void free_arrays();

    void get_view(double &va, double &vb) const;

    void draw_curve(QPainter *painter, double va, double vb,
                    double y_min, double y_max,
                    double (*eval_fn)(double, const Window *),
                    const QColor &color) const;
};

#endif 
