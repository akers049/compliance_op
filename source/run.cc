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

  sprintf(co.ls_static_dir, "../keyword_static");

  co.write_element_file();
  co.write_part_file();

  sprintf(co.run_dir, "../output");
  // see if the directory exists
  struct stat st;
  if (stat(co.run_dir, &st) == -1)
    mkdir(co.run_dir, 0700);

  chdir(co.run_dir);
  getcwd(co.run_dir, MAXLINE);
  chdir("-");

  co.run_ls_dyna(0);
  co.read_displacements();
  co.set_element_stiffness();
  co.compute_sensitivities();
  for(unsigned int i = 0; i < co.N; i ++)
    std::cout << co.sensitivities[i] << std::endl;

  return 0;
}
