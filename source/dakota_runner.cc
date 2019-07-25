#include <iostream>

#include "compliance_opt.h"
#include <stdio.h>
#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>



using namespace std;

int main()
{
  Comp_Op::compliance_opt co;
  getcwd(co.run_dir, MAXLINE);
  co.initialize("../keyword_static");
  unsigned int iter = co.read_dakota_param_input();

  std::cout << "ITERATATION : " << iter << std::endl;
  co.iterate(iter, true);
  co.write_vtk(iter);
  co.set_stiffness_matrix();
  co.compute_hessian();
  co.write_dakota_param_output();

  return 0;
}
