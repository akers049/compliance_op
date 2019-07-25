#include <iostream>
#include <string.h>
#include "compliance_opt.h"
#include <unistd.h>

#define MAXLINE 1024
using namespace std;

int main()
{
  Comp_Op::compliance_opt co;
  char pwd[MAXLINE];
  getcwd(pwd, MAXLINE);
  strcpy(co.run_dir, pwd);
  strcat(co.run_dir, "/workdir_init");

  co.write_dakota_input("keyword_static");
  char command[MAXLINE];
  strcpy(command, "~/dakota/dakota-6.9.0.src/build/src/dakota -input compliance_op_dakota.in");

  system(command);

  return 0;
}
