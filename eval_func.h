#include "Header.h"
static Benchmarks* bench=NULL;
void set_func(int funcID);
void set_data_dir(char *new_data_dir);
double eval_sol(double *x);
void free_func(void);
void next_run(void);
int get_min_region(void);
int get_max_region(void);
