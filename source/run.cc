#include <iostream>

#include "compliance_opt.h"
#include <stdio.h>
#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>



using namespace std;

int main()
{
  char fileName[MAXLINE];
  std::cout << "Please enter an input file: " << std::endl;
  std::cin >> fileName;

  Comp_Op::compliance_opt co;
  co.read_input_file(fileName);
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
/*
//   static
  co.dynamicFlag = false;
  std::fill(co.rho.begin(), co.rho.end(), 0.5);
  co.update_element_rhos();
  for(unsigned int i = 0; i < 20; i ++)
  {
    std::cout << "ITERATATION : " << i << std::endl;
    co.iterate(i);
    co.write_vtk(i);
    std::cout << co.compute_objective() << std::endl;
  }
  exit(-1);
//
 */
  // dynamic

  // do one initialization static iteration
  co.dynamicFlag = false;
  std::cout << "\nINITIAL STATIC ITERATION" << std::endl;
  std::fill(co.rho.begin(), co.rho.end(), 1.0);
  co.update_element_rhos();
  co.iterate(0);

  std::fill(co.rho.begin(), co.rho.end(), 0.5);
  co.update_element_rhos();
  co.setup_dh_drho_vector();

  getcwd(co.dynamic_init_dir, MAXLINE);
  strcat(co.dynamic_init_dir, "/iter_0");

  co.dynamicFlag = true;
  for(unsigned int i = 1; i < 3000; i ++)
  {
    clock_t timer = clock();

    std::cout << "\nITERATATION : " << i << std::endl;
    co.iterate(i);
    std::cout << "OBJECTIVE : " << co.compliance << std::endl;
    std::cout << "Change : " << co.change << std::endl;
    timer = clock() - timer;
    std::cout << "  Time For Iteration : " << ((float)timer)/CLOCKS_PER_SEC << " seconds.\n";



    if(co.change < 0.015)
    {
      std::cout << "\nSTOPPING CRITERIA MET. EXITING." << std::endl;
      break;
    }
  }

  return 0;
}
