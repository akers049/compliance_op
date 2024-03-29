#ifndef MY_COMPLIANCE_OPT_H
#define MY_COMPLIANCE_OPT_H

#include <string>
#include <fstream>

#include <stdlib.h>
#include <math.h>
#include <vector>
#include <cstdlib>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>
#include "MMASolver.h"

#define DIM 2

#define MAXLINE 1024

namespace Comp_Op
{

  void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                        int const maxSize, int* const endOfFileFlag, bool trimFlag = true);

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
     float compute_sensitivity(const std::vector<float> &u);
     float compute_elm_compliance(const std::vector<float> &u);
     double compute_elm_compliance_dynamic(const std::vector<double> &u);
     void add_neighbor( ElmData* nextElm, float dist)
       {neighborList.push_back(std::make_pair(nextElm, dist));}
     float dist(ElmData* otherElm);
     void update_rho(float rho_new, double k1 = 1.0/3.0, double k2 = 3.0);
     void set_vol(float vol){volume_ = vol;};
     void print_info();

     std::vector< std::pair<ElmData*, float> > neighborList;
     std::vector<int> localToGlobal;
     std::vector< std::vector<float> > elmStiffnessMat;
     std::vector< std::vector<float> > nodePos;
     std::vector< float > center;
     std::vector< unsigned int > node_indicies;

     float volume_;
     float rho_;
     float E_;
     double dE_drho_;
     float sensitivity_;



     double penal_;
     unsigned int elmID;
     unsigned int dim_;
     unsigned int num_vertex;
     unsigned int num_dofs;

  private:
     float t;
  };

  class TimeHistory
  {
  public:
    TimeHistory(unsigned int nnodes, unsigned int ndofs, std::vector< std::vector<int> >  &node_ind, double Tmax);
    ~TimeHistory(){clear_lambda_data();};

    void read_data(char *nodoutFile, bool adjointFlag = false, bool loadSymFlag = false);
    void clear_lambda_data()
    {
      for(unsigned int i = 0; i < lambda.size(); i ++)
        free(lambda[i]);

      lambda_time.clear();
      lambda.clear();

    };

    unsigned int N_nodes;
    unsigned int N_dofs;
    std::vector< std::vector<int> > *node_ind_dofs;

    std::vector<double> t;
    std::vector< std::vector<double> > u;
    std::vector< std::vector<double> > v;
    std::vector< std::vector<double> > a;

    std::vector< double* > lambda;
    std::vector< double >  lambda_time;

    double T_max;

    unsigned int numSteps;
  };

  class compliance_opt
  {

  public:
    compliance_opt();
    ~compliance_opt()
    {
      delete timeHistory;
      delete timeHistory_adjoint;
      delete MMA;
    };

    void initialize(char* static_dir);
    void read_input_file(char* fileName);
    float iterate(unsigned int iter, bool dakota = false);
    void postprocess();
    void update_rho();
    void update_rho_MMA();

    void write_element_file();
    void write_part_file();
    void write_mat_file(const unsigned int iter, char *matFileName, bool adjointFlag = false);
    void write_adjoint_load_file(const unsigned int iter, char *loadFileName);
    void set_element_stiffness();
    void run_ls_dyna(const unsigned int iter, bool adjointFlag = false);

    void read_node_data();
    void read_displacements();

    void compute_sensitivities();
    void filter_sensitivities();
    void set_stiffness_matrix();
    void set_mass_matrix();
    void read_mass_matrix();
    void compute_hessian();

    void update_element_vols();
    void update_element_rhos();
    void update_element_neighbors();

    void integrate_dynamic_lambda();
    void setup_dh_drho_vector();
    void get_next_dh_drho(std::vector<std::vector< double > > &dh_drho_j,unsigned int elmID);
    void get_next_dh_drho_time(std::vector<double> &dh_drho_j,unsigned int elmID, unsigned int timeIndex);
    void compute_dynamic_sensitivity_fast();
    void compute_dynamic_dhdrho();
    void integrate_dynamic_sensitivity();
    void compute_dynamic_sensitivity_time_discrete();
    float compute_dynamic_objective();
    float compute_dynamic_disp_obejctive();


    void write_vtk(unsigned int iter);

//    void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
//                          int const maxSize, int* const endOfFileFlag,
//                          bool trimFlag = true);

    float compute_objective();

    void write_dakota_input(char *static_dir);
    unsigned int read_dakota_param_input();
    void write_dakota_param_output();

    MMASolver *MMA;
    TimeHistory *timeHistory;
    TimeHistory *timeHistory_adjoint;

    std::vector<float> rho;
    std::vector<float> old_rho;
    std::vector<double> sensitivities;
    std::vector<double> sensitivities_filtered;

    std::vector<float> u;
    std::vector<std::vector<int> > node_ind_dofs;
    std::vector< ElmData > ElmDatas;
    std::vector<unsigned int > fixedNodes;
    std::vector<std::vector<float> > node_pos;

    Eigen::SparseMatrix<double> K;
    Eigen::SparseMatrix<double> M;
    std::vector<std::vector<double> > hessian;
    std::map< std::pair<unsigned int, unsigned int>, std::vector<double> > dh_drho;
    std::vector<std::vector<std::vector <double> > > dh_drho_vector;
    std::vector< std::vector<unsigned int> > dh_drho_indicies;


    char ls_static_dir[MAXLINE];
    char run_dir[MAXLINE];
    char iter_dir[MAXLINE];
    char dynamic_init_dir[MAXLINE];

    float R_filter;
    float volFrac;
    float rhoMin;
    double penal;

    float density;
    float poisson;

    float V_tot;

    bool firstFlag;
    bool dynamicFlag;
    bool cleanFlag;
    bool loadSymFlag;
    bool MMA_flag;
    bool dispObjectiveFlag;

    unsigned int N;
    unsigned int N_nodes;
    unsigned int N_dofs;

    unsigned int numIntegrationSteps;
    double k1, k2;
    double dispPower;

    float compliance;
    float dispObjective;

    double F_max;
    double T_max;
    std::vector<unsigned int> forcedDofs;
    float change;

    double objectiveFactor;

  private:
    unsigned int internal_iter;

  };


}



#endif
