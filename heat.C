#include <cpgplot.h>
void heat()
{
  int logarithmic=0;
  int itf = 0;
  float contrast=1.0;
  float brightness=0.5;
  if (logarithmic)
    itf = 1;

  cpgsitf (itf);

    float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
    float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
    
    cpgctab (heat_l, heat_r, heat_g, heat_b, 5, contrast, brightness);
}

