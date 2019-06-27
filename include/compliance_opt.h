#ifndef MY_COMPLIANCE_OPT_H
#define MY_COMPLIANCE_OPT_H

#include <string>
#include <fstream>

#include <stdlib.h>
#include <math.h>
#include <vector>
#include <cstdlib>

#define DIM 2

#define MAXLINE 1024

namespace Comp_Op
{
class ElmData
  {

  public:
     ElmData(unsigned int dim, double pen);
     ~ElmData(){};

     void set_elm_stiffness(std::vector<unsigned int> &local_inds,
                       std::vector<std::vector<int> > &node_ind_dofs,
                       std::vector<float> &stiff_values);
     float compute_sensitivity(const std::vector<float> &u, float rho_e);



     std::vector<unsigned int> localToGlobal;
     std::vector< std::vector<float> > elmStiffnessMat;

     double penal_;
     unsigned int dim_;
     unsigned int num_vertex;
     unsigned int num_dofs;
  };

  class compliance_opt
  {

  public:
    compliance_opt();
    ~compliance_opt(){};

    void read_input_file(char* fileName);

    void write_element_file();
    void write_part_file();
    void write_mat_file(const unsigned int iter, char *matFileName);
    void set_element_stiffness();
    void run_ls_dyna(const unsigned int iter);

    void read_node_inds();
    void read_displacements();

    void compute_sensitivities();


    void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                          int const maxSize, int* const endOfFileFlag,
                          bool trimFlag = true);

    std::vector<float> rho;
    std::vector<float> old_rho;
    std::vector<float> sensitivities;
    std::vector<float> u;
    std::vector<std::vector<int> > node_ind_dofs;
    std::vector< ElmData > ElmDatas;

    char ls_static_dir[MAXLINE];
    char run_dir[MAXLINE];
    char iter_dir[MAXLINE];

    float volFrac;
    double rhoMin;
    double penal;

    float density;
    float poisson;

    bool firstFlag;

    unsigned int N;
    unsigned int N_nodes;
    unsigned int N_dofs;

  };


}



#endif
