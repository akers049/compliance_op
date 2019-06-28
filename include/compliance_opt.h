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
     ElmData(unsigned int dim, double pen, unsigned int eID,
         std::vector<unsigned int> &nodeIds,
         std::vector<std::vector<float> > &globalNodePos);
     ~ElmData(){};

     void set_elm_stiffness(std::vector<unsigned int> &local_inds,
                       std::vector<std::vector<int> > &node_ind_dofs,
                       std::vector<float> &stiff_values);
     float compute_sensitivity(const std::vector<float> &u, float rho_e);
     void add_neighbor( ElmData* nextElm, float dist)
       {neighborList.push_back(std::make_pair(nextElm, dist));}
     float dist(ElmData* otherElm);
     void update_rho(float rho_new){rho_ = rho_new;}
     void set_vol(float vol){volume_ = vol;};

     std::vector< std::pair<ElmData*, float> > neighborList;
     std::vector<unsigned int> localToGlobal;
     std::vector< std::vector<float> > elmStiffnessMat;
     std::vector<std::vector <float> > nodePos;
     std::vector< float > center;
     std::vector< unsigned int > node_indicies;

     float volume_;
     float rho_;
     float sensitivity_ = 0.0;
     double penal_;
     unsigned int elmID;
     unsigned int dim_;
     unsigned int num_vertex;
     unsigned int num_dofs;
  };

  class compliance_opt
  {

  public:
    compliance_opt();
    ~compliance_opt(){};

    void initialize(char* static_dir);
    float iterate(unsigned int iter);
    void postprocess();
    void update_rho();
    void read_input_file(char* fileName);

    void write_element_file();
    void write_part_file();
    void write_mat_file(const unsigned int iter, char *matFileName);
    void set_element_stiffness();
    void run_ls_dyna(const unsigned int iter);

    void read_node_data();
    void read_displacements();

    void compute_sensitivities();
    void filter_sensitivities();

    void update_element_vols();
    void update_element_rhos();
    void update_element_neighbors();
    void write_vtk(unsigned int iter);

    void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                          int const maxSize, int* const endOfFileFlag,
                          bool trimFlag = true);

    float compute_objective(){return u[node_ind_dofs[220][1]]; }

    std::vector<float> rho;
    std::vector<float> old_rho;
    std::vector<float> sensitivities;
    std::vector<float> sensitivities_filtered;

    std::vector<float> u;
    std::vector<std::vector<int> > node_ind_dofs;
    std::vector< ElmData > ElmDatas;
    std::vector<unsigned int > fixedNodes;
    std::vector<std::vector<float> > node_pos;

    char ls_static_dir[MAXLINE];
    char run_dir[MAXLINE];
    char iter_dir[MAXLINE];

    float R_filter;
    float volFrac;
    float rhoMin;
    double penal;

    float density;
    float poisson;

    float V_tot;

    bool firstFlag;

    unsigned int N;
    unsigned int N_nodes;
    unsigned int N_dofs;

  private:
    unsigned int internal_iter = 0;

  };


}



#endif
