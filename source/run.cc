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

  co.initialize("../keyword_static");

  sprintf(co.run_dir, "../output");
  // see if the directory exists
  struct stat st;
  if (stat(co.run_dir, &st) == -1)
    mkdir(co.run_dir, 0700);

  chdir(co.run_dir);
  getcwd(co.run_dir, MAXLINE);
  chdir("-");

  co.iterate(0);


  return 0;
}
