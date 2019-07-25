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
  co.initialize("keyword_static");

//  sprintf(co.run_dir, "../output");
//  // see if the directory exists
//  struct stat st;
//  if (stat(co.run_dir, &st) == -1)
//    mkdir(co.run_dir, 0700);
//
//  chdir(co.run_dir);
//  getcwd(co.run_dir, MAXLINE);
//  chdir("-");

  // static
//  for(unsigned int i = 0; i < 20; i ++)
//  {
//    std::cout << "ITERATATION : " << i << std::endl;
//    co.iterate(i);
//    co.write_vtk(i);
//    std::cout << co.compute_objective() << std::endl;
//
//
  // dynamic

  // do one initialization static iteration
  co.dynamicFlag = false;
  std::cout << "INITIAL STATIC ITERATION" << std::endl;
  std::fill(co.rho.begin(), co.rho.end(), 1.0);
  co.update_element_rhos();
  co.iterate(0);

  std::fill(co.rho.begin(), co.rho.end(), 0.5);
  co.update_element_rhos();

  getcwd(co.dynamic_init_dir, MAXLINE);
  strcat(co.dynamic_init_dir, "/iter_0");

  co.dynamicFlag = true;
  for(unsigned int i = 1; i < 20; i ++)
  {
    std::cout << "ITERATATION : " << i << std::endl;
    co.iterate(i);
    std::cout << co.compliance << std::endl;

  }

  return 0;
}
