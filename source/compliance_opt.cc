#ifndef MY_COMPLIANCE_OPT_C
#define MY_COMPLIANCE_OPT_C

#include "compliance_opt.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <stdlib.h>
#include <cmath>
#include <unistd.h>

#include <gsl/gsl_poly.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <time.h>

double interpolate( std::vector<double> &xData, std::vector<double> &yData, double x, bool extrapolate )
{
   int size = xData.size();

   int i = 0;                                                                  // find left end of interval for interpolation
   if ( x >= xData[size - 2] )                                                 // special case: beyond right end
   {
      i = size - 2;
   }
   else
   {
      while ( x > xData[i+1] ) i++;
   }
   double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
   if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
   {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
   }

   double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

   return yL + dydx * ( x - xL );                                              // linear interpolation
}

double loadingFunction(double t, double F_max, double dT)
{
  double f = ((F_max*2.0/dT)*fabs(-t - dT/2.0) - F_max);

  return f;
}

double dloadingFunction(double t, double F_max, double dT)
{
  double f = 0.0;
  if(-t < dT/2.0)
    f = -F_max*2.0/dT;
  else
    f = F_max*2.0/dT;

  return f;
}


// gsl integration functions
int
func (double t, const double y[], double f[],
      void *params)
{

  // unpack the stiffness and mass matrix
  std::vector< void * > * PD =
                       (std::vector< void * > *) params;
  Eigen::SparseMatrix<double> * K = (Eigen::SparseMatrix<double> *) (*PD)[0];
  Eigen::SparseMatrix<double> * M = (Eigen::SparseMatrix<double> *) (*PD)[1];
  double F_max = *(double *) (*PD)[2];
  double T_max = *(double *) (*PD)[3];
  std::vector<unsigned int> *loadDofs = (std::vector<unsigned int> *) (*PD)[4];

  unsigned int N_dofs = K->rows();
  // zero out f
  for(unsigned int i = 0; i < 2.0*N_dofs; i ++)
    f[i] = 0.0;

  // Add in the stiffness matrix contributions
  for (unsigned int i = 0; i < K->outerSize(); ++i)
    for (Eigen::SparseMatrix<double>::InnerIterator it((*K),i); it; ++it)
    {
      f[it.row()] += it.value()*y[it.col() + N_dofs];
    }

  // Add in the mass matrix contributions
  for (unsigned int i = 0; i < M->outerSize(); ++i)
    for (Eigen::SparseMatrix<double>::InnerIterator it((*M),i); it; ++it)
    {
      f[it.row() + N_dofs] += -it.value()*y[it.col()];
    }

  // add in the forcing contribution
  double f_value = loadingFunction(t, F_max, T_max);
  for(unsigned int i = 0 ; i < loadDofs->size(); i ++)
    f[(*loadDofs)[i]] += f_value;

  // Scale by the mass
  for (unsigned int i = 0; i < M->outerSize(); ++i)
    for (Eigen::SparseMatrix<double>::InnerIterator it((*M),i); it; ++it)
    {
      double nextVal = 1.0/it.value();
      f[it.row()] *= nextVal;
      f[it.row() + N_dofs] *= nextVal;
    }



  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy,
     double dfdt[], void *params)
{
  // unpack the stiffness and mass matrix
  std::vector< void * > * PD =
                       (std::vector< void * > *) params;
  Eigen::SparseMatrix<double> * K = (Eigen::SparseMatrix<double> *) (*PD)[0];
  Eigen::SparseMatrix<double> * M = (Eigen::SparseMatrix<double> *) (*PD)[1];
  double F_max = *(double *) (*PD)[2];
  double T_max = *(double *) (*PD)[3];
  std::vector<unsigned int> *loadDofs = (std::vector<unsigned int> *) (*PD)[4];

  unsigned int N_dofs = K->rows();

  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 2*N_dofs, 2*N_dofs);
  gsl_matrix * m = &dfdy_mat.matrix;

  // zero out all the values for now
  for(unsigned int i = 0; i < 2*N_dofs; i ++)
    for(unsigned int j = 0; j < 2*N_dofs; j ++)
      gsl_matrix_set(m, i, j, 0.0);

  // Add in the stiffness matrix contributions
  for (unsigned int i = 0; i < K->outerSize(); ++i)
    for (Eigen::SparseMatrix<double>::InnerIterator it((*K),i); it; ++it)
    {
      gsl_matrix_set(m, it.row(), it.col() + N_dofs, it.value()/M->coeffRef(it.row(), it.row()));
    }

  // Add in the mass matrix contributions
  for (unsigned int i = 0; i < M->outerSize(); ++i)
    for (Eigen::SparseMatrix<double>::InnerIterator it((*M),i); it; ++it)
    {
      gsl_matrix_set(m, it.row() + N_dofs, it.col() + N_dofs, -1.0);
    }

  // zero out dfdt
  for(unsigned int i = 0; i < 2.0*N_dofs; i ++)
    dfdt[i] = 0.0;

  // add in the forcing contribution to dfdt
  double df_value = dloadingFunction(t, F_max, T_max);
  for(unsigned int i = 0 ; i < loadDofs->size(); i ++)
    dfdt[(*loadDofs)[i]] += df_value;


  // Scale by the mass
  for (unsigned int i = 0; i < M->outerSize(); ++i)
    for (Eigen::SparseMatrix<double>::InnerIterator it((*M),i); it; ++it)
    {
      double nextVal = 1.0/it.value();
      dfdt[it.row()] *= nextVal;
    }

  return GSL_SUCCESS;
}


namespace Comp_Op
{


/******************** time_history member functions **********************/
  TimeHistory::TimeHistory(unsigned int nnodes, unsigned int ndofs, std::vector< std::vector<int> >  &node_ind, double Tmax)
  {
    N_dofs = ndofs;
    N_nodes = nnodes;
    node_ind_dofs = &(node_ind);
    numSteps = 0;
    T_max = Tmax;
  }

  void TimeHistory::read_data(char *nodoutFile, bool adjointFlag, bool loadSymFlag)
  {
    double maxVal = 0.0;
    numSteps = 0;
    u.clear();
    v.clear();
    a.clear();
    t.clear();

    FILE* dispFile;
    dispFile = fopen(nodoutFile, "r");
    if(dispFile == NULL)
    {
      std::cout << "Error opening displacement file, Exiting.\n" << std::endl;
      exit(1);
    }

    clock_t timer = clock();

    int dummy;
    std::vector<int> tmp(DIM);
    char nextLine[MAXLINE];
    unsigned int count = 0;
    char *nextDataString;
    getNextDataLine(dispFile, nextLine, MAXLINE, &dummy, false);
    while(dummy == 0)
    {
      if(strncmp(nextLine, " n o d a l", 10) == 0)
      {
        count ++;
        if(count%2 == 1)
        {
          numSteps ++;
          // these are the displacements, not the rotations
          u.push_back(std::vector<double> (N_dofs, 0.0));
          v.push_back(std::vector<double> (N_dofs, 0.0));
          a.push_back(std::vector<double> (N_dofs, 0.0));
          t.push_back(0.0);

          std::vector<double> nextVals(9);
          unsigned int nextNodeNumber;

          // save time
          char *timeChar = &(nextLine[strlen(nextLine) - 16]);
          unsigned int valuesWritten = sscanf(timeChar, "%lf", &(t[numSteps-1]));
          if(valuesWritten != 1)
          {
            std::cout << "Error reading the time history data. Exiting.\n";
            exit(-1);
          }

          getNextDataLine(dispFile, nextLine, MAXLINE, &dummy, true); //blank lines
          for(unsigned int i = 0; i < N_nodes; i ++)
          {
            getNextDataLine(dispFile, nextLine, MAXLINE, &dummy, true); //blank lines

            // read data and node number
//            valuesWritten = sscanf(nextLine, "%u %lf %lf %lf %lf %lf %lf %lf %lf %lf",
//                &nextNodeNumber, &(nextVals[0]), &(nextVals[1]), &(nextVals[2]),
//                                 &(nextVals[3]), &(nextVals[4]), &(nextVals[5]),
//                                 &(nextVals[6]), &(nextVals[7]), &(nextVals[8]));
//            if(valuesWritten != 10)
//            {
//              std::cout << "Error reading the time history data. Exiting.\n";
//              exit(-1);
//            }

            valuesWritten = sscanf(nextLine, "%u", &nextNodeNumber);
            if(valuesWritten != 1)
            {
              std::cout << "Error reading the time history data. Exiting.\n";
              exit(-1);
            }

            nextDataString = &nextLine[10];
            for(unsigned int j = 0; j < 8; j ++)
            {
              if(nextDataString[8] == 'E')
              {
                valuesWritten = sscanf(nextDataString, "%lf", &(nextVals[j]));
                if(valuesWritten != 1)
                {
                  std::cout << "Error reading the time history data. Exiting.\n";
                  exit(-1);
                }
              }
              else
                nextVals[j] = 0.0;

              nextDataString = &(nextDataString[12]);
            }
            if(nextDataString[8] == 'E')
            {
              valuesWritten = sscanf(nextDataString, "%lf", &(nextVals[8]));
              if(valuesWritten != 1)
              {
                std::cout << "Error reading the time history data. Exiting.\n";
                exit(-1);
              }
            }
            else
              nextVals[8] = 0.0;


            for(unsigned int j = 0; j < 9; j ++)
            {
              if(fabs(nextVals[j]) < 1e-15)
                nextVals[j] = 0.0;

              if(numSteps == 3)
              {
                if(fabs(nextVals[j]) > maxVal)
                  maxVal = fabs(nextVals[j]);
              }

            }
            for(unsigned int j = 0; j < DIM; j ++)
            {
              if((*node_ind_dofs)[nextNodeNumber -1][j] == -1)
                continue;

              // put in the data ya know for u v a
              u[numSteps -1][(*node_ind_dofs)[nextNodeNumber -1][j]] = nextVals[j];
              v[numSteps -1][(*node_ind_dofs)[nextNodeNumber -1][j]] = nextVals[j + 3];
              a[numSteps -1][(*node_ind_dofs)[nextNodeNumber -1][j]] = nextVals[j + 6];

            }

          }
        }
      }

      getNextDataLine(dispFile, nextLine, MAXLINE, &dummy, false);
    }
    if(dummy == 1 && numSteps == 0)
    {
      std::cout << "Error reading displacement file. Exiting. \n" << std::endl;
      exit(-1);
    }

    timer = clock() - timer;


    if(adjointFlag)
    {
      double multFactor = (loadSymFlag ? 1.0 : -1.0);
      clear_lambda_data();
      lambda.resize(numSteps);
      lambda_time.resize(numSteps);
      for(unsigned int i = 0; i < numSteps; i++)
      {
        lambda[i] = (double *) calloc(2*N_dofs, sizeof(double));
        lambda_time[i] = T_max - t[numSteps - i - 1];

        for(unsigned int k = 0; k < N_dofs; k++)
        {
          lambda[i][k] = multFactor*v[numSteps - i - 1][k];
          lambda[i][k + N_dofs] = -multFactor*u[numSteps - i - 1][k];
        }
      }
    }

    fclose(dispFile);
    std::cout << "  Sucessfully read in time history file with " << numSteps << " timesteps.\n";
    std::cout << "    Elapsed time to read history file : " << ((float)timer)/CLOCKS_PER_SEC << " seconds.\n";
    std::cout << "    MAX VAL : " << maxVal << std::endl;

  }



/******************** ElmData member functions **********************/

  ElmData::ElmData(unsigned int dim, double pen, unsigned int eID,
      std::vector<unsigned int> &nodeIds,
      std::vector<std::vector<float> > &globalNodePos)
  :
  dim_(dim),
  penal_(pen),
  elmID(eID),
  rho_(0.0)
  {
    num_vertex = pow(2, dim);
    num_dofs = num_vertex*dim;

    elmStiffnessMat.resize(num_dofs, std::vector<float> (num_dofs, 0.0));
    localToGlobal.resize(num_dofs, -1);
    center.resize(DIM, 0.0);
    float multFactor = 1.0/(1.0*num_vertex);
    for(unsigned int i = 0; i < num_vertex; i ++)
    {
      nodePos.push_back(globalNodePos[nodeIds[i]-1]);

      for(unsigned int k = 0; k < DIM; k ++)
        center[k] += multFactor*globalNodePos[nodeIds[i]-1][k];
    }

    volume_ = 1.0;
    node_indicies = nodeIds;
    sensitivity_ = 0.0;
    t = 0.0;
    E_ = 0.0;
    dE_drho_ = 0.0;
  }

  void ElmData::update_rho(float rho_new, double k1, double k2)
  {
    rho_ = rho_new;

    double a = -3.0*(1.0 - k2)/(k1 - k2);
    double b = -a;
    double c = -rho_;

    double roots[3];
    int numRoots = gsl_poly_solve_cubic(a, b, c, &roots[0], &roots[1], &roots[2]);
    for(unsigned int i = 0; i < numRoots; i++)
      if(roots[i] >= 0.0 && roots[i] <= 1.0)
      {
        t = roots[i];
        break;
      }
    double t_sq = t*t;
    E_ = k1*(a*t_sq + b*t) + t_sq*t;

    dE_drho_ = (k1*b*(1.0 - 2.0*t) + 3.0*t_sq)/(b*(1.0 - 2.0*t) + 3.0*t_sq);

  }

  void ElmData::set_elm_stiffness(std::vector<unsigned int> &local_inds,
                                  std::vector<std::vector<int> > &node_ind_dofs,
                                  std::vector<float> &stiff_values)
  {
    elmStiffnessMat.resize(num_dofs, std::vector<float> (num_dofs, 0.0));
    localToGlobal.resize(num_dofs, -1);

    unsigned int count = 0;
    for(unsigned int i = 0; i < num_vertex; i ++)
      for(unsigned int k = 0; k < dim_; k ++)
      {
        localToGlobal[count] = node_ind_dofs[local_inds[i] - 1][k];
        count ++;
      }

    count = 0;
    for(unsigned int i = 0; i < num_dofs; i ++)
      for(unsigned int k = 0; k <= i; k ++)
      {
        elmStiffnessMat[k][i] = stiff_values[count];

        if(i != k)
          elmStiffnessMat[i][k] = elmStiffnessMat[k][i];

        count ++;
      }
  }

  float ElmData::dist(ElmData *otherElm)
  {
    float distance = 0.0;
    for(unsigned int i = 0; i < dim_; i ++)
      distance += (center[i]- otherElm->center[i])*(center[i]- otherElm->center[i]);

    distance = sqrt(distance);

    return distance;
  }

  float ElmData::compute_sensitivity(const std::vector<float> &u)
  {
    sensitivity_ = 0.0;
    for(unsigned int i = 0; i < num_dofs; i ++)
      for(unsigned int k = 0; k < num_dofs; k ++)
      {
        if(localToGlobal[i] < 0 || localToGlobal[k] < 0)
          continue;
        else
          sensitivity_ += u[localToGlobal[i]]*elmStiffnessMat[i][k]*u[localToGlobal[k]];
      }

    sensitivity_ *= (-penal_/rho_);

    return  sensitivity_;
  }

  double ElmData::compute_elm_compliance_dynamic(const std::vector<double> &u)
  {
    double elm_compliance_dub = 0.0;
    for(unsigned int i = 0; i < num_dofs; i ++)
      for(unsigned int k = 0; k < num_dofs; k ++)
      {
        if(localToGlobal[i] < 0 || localToGlobal[k] < 0)
          continue;
        else
          elm_compliance_dub += E_*u[localToGlobal[i]]*elmStiffnessMat[i][k]*u[localToGlobal[k]];
      }

    return elm_compliance_dub;
  }


  float ElmData::compute_elm_compliance(const std::vector<float> &u)
  {
    float elm_compliance = 0.0;
    for(unsigned int i = 0; i < num_dofs; i ++)
      for(unsigned int k = 0; k < num_dofs; k ++)
      {
        if(localToGlobal[i] < 0 || localToGlobal[k] < 0)
          continue;
        else
          elm_compliance += u[localToGlobal[i]]*elmStiffnessMat[i][k]*u[localToGlobal[k]];
      }

    return elm_compliance;
  }

  void ElmData::print_info()
  {
    std::cout << "     ELMID : " << elmID << std::endl;
    std::cout << "     CENTER : " << center[0] << " " << center[1] << std::endl;
    std::cout << "     SENSITIVITY : " << sensitivity_ << std::endl;
    std::cout << "     DISP LOCAL : " << std::endl;
    for(unsigned int i = 0; i < num_dofs; i++)
      std::cout << "    " << localToGlobal[i] << std::endl;
    std::cout << "     STIFFNESS MAT LOCAL : " << std::endl;
    for(unsigned int i = 0; i < num_dofs; i++)
    {
      for(unsigned int k = 0; k < num_dofs; k++)
        std::cout << elmStiffnessMat[i][k] << " ";
      std::cout << std::endl;
    }
    std::cout << "\n";

  }

  /******************** compliance_opt member functions **********************/

  compliance_opt::compliance_opt()
  {
    volFrac = 0.0;
    rhoMin = 0.0;
    penal = 0.0;

    density = 0;
    poisson = 0;
    N = 0;
    N_nodes = 0;
    N_dofs = 0;
    R_filter = 0.0;

    V_tot = 0.0;

    firstFlag = true;


    compliance = 0.0;

    T_max = 0.0;
    F_max = 0.0;
    MMA = NULL;
    timeHistory = NULL;
    timeHistory_adjoint = NULL;

    internal_iter = 0;
    dynamicFlag = true;

    k1 = 1.0/3.0;
    k2 = 3.0;
    change = 10000;

    numIntegrationSteps = 100;
    cleanFlag = true;
    loadSymFlag = true;
    MMA_flag = true;

    dispPower = 4.0;

  }

  void compliance_opt::read_input_file(char* fileName)
  {
    FILE* fid;
    int endOfFileFlag;
    char nextLine[MAXLINE];

    int valuesWritten;
    bool fileReadErrorFlag = false;
    unsigned int dummy;


    fid = std::fopen(fileName, "r");
    if (fid == NULL)
    {
      std::cout << "Unable to open file \"" << fileName  << "\"" <<  std::endl;
      fileReadErrorFlag = true;
    }
    else
    {

      // Read in the number of elements, number of nodes, and number of dofs
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u %u %u", &N, &N_nodes, &N_dofs);
      if (valuesWritten != 3)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the density and poisson ratio
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%f %f", &density, &poisson);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in constrained volume fraction, the minimum rho, the filter
      // radius, and the penalization factor.
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%f %f %f %lf", &volFrac, &rhoMin, &R_filter, &penal);
      if(valuesWritten != 4)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      //reading forcing info: T_max and F_max
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lf %lf", &T_max, &F_max);
      if(valuesWritten != 2)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }

      // read in the forced dof
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u", &dummy);
      if(valuesWritten != 1)
      {
        fileReadErrorFlag = true;
        goto fileClose;
      }
      forcedDofs.push_back(dummy);

      // read in k1 and k2 for dynamics
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%lf %lf", &k1, &k2);
      if(valuesWritten != 2)
      {
        std::cout << "k1 and k2 not set. Using default values of 1/3 and 3, respectively. \n";
      }

      // read in the number of integration steps for dynamics
      getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
      valuesWritten = sscanf(nextLine, "%u", &numIntegrationSteps);
      if(valuesWritten != 1)
      {
        std::cout << "Number of integration steps not set. Using default value of 100.\n";
      }


      fileClose:
      {
        fclose(fid);
      }
    }

    if (fileReadErrorFlag)
    {
      // default parameter values
      std::cout << "  Error reading input file, Exiting.\n" << std::endl;
      exit(1);
    }
    else
      std::cout << "  Input file successfully read" << std::endl;
  }

  void compliance_opt::initialize(char* static_dir)
  {

    node_ind_dofs.resize(N_nodes, std::vector<int>(DIM));
    node_pos.resize(N_nodes, std::vector<float>(DIM));

    u.resize(N_dofs);
    sensitivities.resize(N);
    sensitivities_filtered.resize(N, 0.0);

    rho.resize(N);
    std::fill(rho.begin(), rho.end(), volFrac);
    K.resize(N_dofs, N_dofs);
    M.resize(N_dofs, N_dofs);

    char currentDir[MAXLINE];
     getcwd(currentDir, MAXLINE);

    chdir(static_dir);
    getcwd(ls_static_dir, MAXLINE);


    read_node_data();

    timeHistory = new TimeHistory(N_nodes, N_dofs, node_ind_dofs, T_max);
    timeHistory_adjoint = new TimeHistory(N_nodes, N_dofs, node_ind_dofs, T_max);

    MMA = new MMASolver(N, 1);

    write_element_file();

    update_element_rhos();
    write_part_file();

    update_element_neighbors();

  }

  void compliance_opt::postprocess()
  {
    // just gotta do this once
    if(internal_iter == 0)
    {
      update_element_vols();
    }

    if(!dynamicFlag)
    {
      set_element_stiffness();
      read_displacements();
      compute_sensitivities();
      compute_objective();
    }
    else
    {
//      set_stiffness_matrix();
//      set_mass_matrix();

      char dispFileName[MAXLINE];
      strcpy(dispFileName, iter_dir);
      strcat(dispFileName, "/nodout");

      char dispFileName_2[MAXLINE];
      strcpy(dispFileName_2, dispFileName);
      strcat(dispFileName_2, "_forward");
      timeHistory->read_data(dispFileName);

      // adjoint run
      if(!loadSymFlag)
        run_ls_dyna(0, true);

      timeHistory_adjoint->read_data(dispFileName, true, loadSymFlag);

//      integrate_dynamic_lambda();
      compute_dynamic_dhdrho();
      integrate_dynamic_sensitivity();
      compute_dynamic_objective();
    }

    filter_sensitivities();

  }

  float compliance_opt::iterate(unsigned int iter, bool dakota)
  {
    char nextDir[MAXLINE];

    if(dakota == false)
    {
      char iter_char[8];
      sprintf(iter_char, "%u", iter);
      strcpy(nextDir, "iter_");
      strcat(nextDir, iter_char);

      struct stat st;
      if (stat(nextDir, &st) == -1)
        mkdir(nextDir, 0700);

      chdir(nextDir);
    }

    if((internal_iter != 0 && !dynamicFlag) ||
         (dynamicFlag && internal_iter > 1))
    {
      if(MMA_flag)
        update_rho_MMA();
      else
        update_rho();

      update_element_rhos();
    }

    run_ls_dyna(iter);
    postprocess();
    write_vtk(iter);

    if(dakota == false)
    {
      chdir("..");

      if(cleanFlag)
      {
        char rmCommand[MAXLINE];
        strcpy(rmCommand, "rm -r ");
        strcat(rmCommand, nextDir);
        system(rmCommand);
      }
    }
    internal_iter++;

    return compliance;

  }

  void compliance_opt::run_ls_dyna(const unsigned int iter, bool adjointFlag)
  {
//    chdir(run_dir);
    char matFileName[MAXLINE];
    write_mat_file(iter, matFileName, adjointFlag);

    char command[MAXLINE];

//    sprintf(command, "cd %s;", iter_dir);
    strcpy(command, "module load bhatta/ls-dyna-double-smp/r1010;");

    char runCommand[MAXLINE];
    sprintf(runCommand,
        "/groups/bhatta/software/ls-dyna_smp_d_r1010_x64_redhat5_ifort160/ls-dyna_smp_d_r1010_x64_redhat5_ifort160 i=%s >>ResultOutput.txt",
        matFileName);

    strcat(command, runCommand);

//    if(adjointFlag)
//    {
//      char mv_command[MAXLINE];
//      strcpy(mv_command, "mv nodout nodout_forward;");
//      system(mv_command);
//    }

    system(command);

  }

  void compliance_opt::read_node_data()
    {
      // first read the ls_info.k file to get the constrained nodes.
      char infoFileName[MAXLINE];
      strcpy(infoFileName, ls_static_dir);
      strcat(infoFileName, "/ls_info.k");

      FILE* infoFile;
      infoFile = fopen(infoFileName, "r");
      if(infoFile == NULL)
      {
        std::cout << "Error opening node file " << infoFileName << ", Exiting.\n" << std::endl;
        exit(1);
      }

      // read in the fixed nodes
      unsigned int valuesWritten;
      char nextLine[MAXLINE];
      int fileEndFlag;
      getNextDataLine(infoFile, nextLine, MAXLINE, &fileEndFlag, false);
      while(fileEndFlag == false)
      {
        if(strncmp(nextLine, "*SET_NODE_LIST_TITLE", 20) == 0)
        {
          // skip this one
          getNextDataLine(infoFile, nextLine, MAXLINE, &fileEndFlag, false);

          // line with sid
          unsigned int sid;
          getNextDataLine(infoFile, nextLine, MAXLINE, &fileEndFlag, false);
          if(fileEndFlag == true)
          {
            std::cout << "Error 1 reading ls_info.k file, Exiting.\n" << std::endl;
            exit(-1);
          }
          valuesWritten = sscanf(nextLine, "%u", &sid);
          if(valuesWritten != 1)
          {
            std::cout << nextLine << std::endl;
            std::cout << "Error 2 reading ls_info.k file, Exiting.\n" << std::endl;
            exit(-1);
          }
          if(sid == 1)
            break;
        }

        getNextDataLine(infoFile, nextLine, MAXLINE, &fileEndFlag, false);
      }
      if(fileEndFlag == true)
      {
        std::cout << "Error 3 reading ls_info.k file, Exiting.\n" << std::endl;
        exit(1);
      }

      std::vector<unsigned int> node_nums(8);
      getNextDataLine(infoFile, nextLine, MAXLINE, &fileEndFlag, false);
      while(fileEndFlag == false)
      {
        valuesWritten = sscanf(nextLine,
            "%u %u %u %u %u %u %u %u",
            &node_nums[0], &node_nums[1], &node_nums[2], &node_nums[3],
            &node_nums[4], &node_nums[5], &node_nums[6], &node_nums[7]);

        if(valuesWritten == 0)
          break;
        else
        {
          for(unsigned int i = 0 ; i < valuesWritten; i ++)
          {
            if(node_nums[i] == 0)
              continue;

            fixedNodes.push_back(node_nums[i] - 1);
          }
        }

        getNextDataLine(infoFile, nextLine, MAXLINE, &fileEndFlag, false);
      }
      if(fileEndFlag == true)
      {
        std::cout << "Error 4 reading ls_info.k file, Exiting.\n" << std::endl;
        exit(1);
      }

      unsigned int count = 0;
      for(unsigned int i = 0; i < N_nodes; i ++)
      {
        bool found;
        std::vector<unsigned int>::iterator it = std::find(fixedNodes.begin(), fixedNodes.end(), i);
        if(it == fixedNodes.end())
          found = false;
        else
          found = true;

        for(unsigned int k = 0; k < DIM; k ++)
        {
          if(found == true)
            node_ind_dofs[i][k] = -1;
          else
          {
            node_ind_dofs[i][k] = count;
            count ++;
          }
        }
      }

      getNextDataLine(infoFile, nextLine, MAXLINE, &fileEndFlag, false);
      while(fileEndFlag == false)
      {
        if(strncmp(nextLine, "*NODE", 5) == 0)
          break;

        getNextDataLine(infoFile, nextLine, MAXLINE, &fileEndFlag, false);
      }
      if(fileEndFlag == true)
      {
        std::cout << "Error 5 reading ls_info.k file, Exiting.\n" << std::endl;
        exit(1);
      }

      char* tokPtr;
      unsigned int nid;
      for(unsigned int i = 0; i < N_nodes; i ++)
      {
        getNextDataLine(infoFile, nextLine, MAXLINE, &fileEndFlag, false);
        if(fileEndFlag == true)
        {
          std::cout << "Error 6 reading ls_info.k file, Exiting.\n" << std::endl;
          exit(1);
        }
        tokPtr = strtok(nextLine, " ");
        valuesWritten = sscanf(tokPtr, "%u", &nid);
        if(valuesWritten != 1)
        {
          std::cout << "Error 7 reading ls_info.k file, Exiting.\n" << std::endl;
          exit(1);
        }

        for(unsigned int k = 0; k < DIM; k ++)
        {
          tokPtr = strtok(NULL, " ");
          valuesWritten = sscanf(tokPtr, "%f", &node_pos[nid - 1][k]);
          if(valuesWritten != 1)
          {
            std::cout << "Error 8 reading ls_info.k file, Exiting.\n" << std::endl;
            exit(1);
          }
        }
      }
      fclose(infoFile);

      std::cout << "  Sucessfully read node file.\n";


    }

  void compliance_opt::write_element_file()
  {
    unsigned int eid, pid;
    std::vector<unsigned int > n(8);

    chdir(ls_static_dir);

    bool fileReadErrorFlag = false;
    char filename[] = "ls_elm_file_raw.k";
    FILE* fid;
    fid = std::fopen(filename, "r");
    if (fid == NULL)
    {
      std::cout << "Unable to open file \"" << filename  << "\"" <<  std::endl;
      return;
    }

    char filename_write[] = "ls_elm_file.k";
    FILE* fid_w;
    fid_w = std::fopen( filename_write, "w");
    if (fid_w == NULL)
    {
      std::cout << "Unable to open file \"" << filename_write  << "\"" <<  std::endl;
      return;
    }

    fprintf(fid_w, "*KEYWORD\n");
    fprintf(fid_w, "*ELEMENT_SHELL\n");
    fprintf(fid_w, "$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8\n");

    char nextLine[MAXLINE];
    int fileEndFlag = false;
    int valuesWritten = 0;

    getNextDataLine(fid, nextLine,MAXLINE, &fileEndFlag, true);
    getNextDataLine(fid, nextLine,MAXLINE, &fileEndFlag, true);
    valuesWritten = sscanf(nextLine, "%u %u %u %u %u %u %u %u %u %u",
            &eid, &pid, &n[0], &n[1], &n[2], &n[3], &n[4], &n[5], &n[6], &n[7]);
    if (valuesWritten != 10)
    {
      fileReadErrorFlag = true;
    }
    unsigned int count = 0;
    while(fileEndFlag == false && fileReadErrorFlag == false)
    {
      if(count != eid - 1)
      {
        std::cout << "Element with eid :  " << eid << " not found in " << count << " position." << std::endl;
      }
      ElmDatas.push_back(ElmData(DIM, penal, eid, n, node_pos));
      fprintf(fid_w, "%8u%8u%8u%8u%8u%8u%8u%8u%8u%8u\n",
        eid, eid, n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7]);

      getNextDataLine(fid, nextLine,MAXLINE, &fileEndFlag, true);
      valuesWritten = sscanf(nextLine, "%u %u %u %u %u %u %u %u %u %u",
          &eid, &pid, &n[0], &n[1], &n[2], &n[3], &n[4], &n[5], &n[6], &n[7]);
      if (valuesWritten != 10)
      {
        fileReadErrorFlag = true;
        break;
      }
      count ++;
    }
    fprintf(fid_w, "*END");

    fclose(fid);
    fclose(fid_w);

    if (fileReadErrorFlag)
    {
      std::cout << "  Error reading file, Exiting.\n" << std::endl;
      exit(1);
    }
    else
     std::cout << "  Successfully wrote element file" << std::endl;


    chdir(run_dir);

  }

  void compliance_opt::write_part_file()
  {
    chdir(ls_static_dir);

    char filename_write[] = "ls_part_file.k";
    FILE* fid_w;
    fid_w = std::fopen( filename_write, "w");
    if (fid_w == NULL)
    {
      std::cout << "Unable to open file \"" << filename_write  << "\"" <<  std::endl;
      return;
    }
    fprintf(fid_w, "*KEYWORD\n");

    for(unsigned int i = 0; i < N; i ++)
    {
      fprintf(fid_w, "*PART\nplate%u\n", i+1);

      fprintf(fid_w, "%10u%10u%10u         0         1\n",
                       i+1, 1, i+1);
    }

    fprintf(fid_w, "*DATABASE_HISTORY_NODE*\n");
    for(unsigned int i = 1; i <= N_nodes; i ++)
    {
      fprintf(fid_w, "%10u", i);
      if(i%8 == 0)
        fprintf(fid_w, "\n");
    }
    if(N_nodes%8 != 0)
      fprintf(fid_w, "\n");

    fprintf(fid_w, "*END");
    fclose(fid_w);

    std::cout << "  Part File Sucessfully Written\n";

    chdir(run_dir);
  }

  void compliance_opt::write_mat_file(const unsigned int iter, char *matFileName, bool adjointFlag)
  {
    char iter_char[8];
    sprintf(iter_char, "%u", iter);
    strcpy(iter_dir, ".");

    strcpy(matFileName, "ls_mat_");
    strcat(matFileName, iter_char);
    if(adjointFlag)
      strcat(matFileName, "_adjoint");
    strcat(matFileName, ".k");

    FILE* matFile;
    matFile = fopen(matFileName, "w");
    if(matFile == NULL)
    {
      std::cout << "Error opening file, Exiting.\n" << std::endl;
      exit(1);
    }

    fprintf(matFile, "*keyword\n");

    float dummy = 0.0;
    unsigned int dummy_int = 0;
    for(unsigned int i = 0; i < N; i ++)
    {
      fprintf(matFile, "*MAT_ELASTIC_TITLE\nelastic%u\n", i);
      if(dynamicFlag)
        fprintf(matFile, "%10u %9f %9f %9f %9g %9g %9u\n", ElmDatas[i].elmID , (ElmDatas[i].rho_)*density, ElmDatas[i].E_, poisson,
                  dummy, dummy, dummy_int);
      else
        fprintf(matFile, "%10u %9f %9f %9f %9g %9g %9u\n", ElmDatas[i].elmID , 1.0, (float) pow(ElmDatas[i].rho_, penal), poisson,
                  dummy, dummy, dummy_int);

    }
    fprintf(matFile, "*include_path_relative\n../keyword_static\n");
    fprintf(matFile, "*include\n");
    char ls_info_path[MAXLINE];

    if(dynamicFlag == false)
      sprintf(ls_info_path, "ls_info.k");
    else if(adjointFlag == false)
      sprintf(ls_info_path, "ls_info_dynamic.k");
    else if(adjointFlag == true)
      sprintf(ls_info_path, "ls_info_dynamic_adjoint.k");

    unsigned int lines = std::floor((1.0*strlen(ls_info_path))/78.0);
    for(unsigned int i = 0; i < lines; i ++)
    {
      char nextLine[78];
      strncpy(nextLine, &(ls_info_path[i*78]), 78);
      fprintf(matFile, "%s +\n", nextLine);
    }
    fprintf(matFile, "%s\n", &(ls_info_path[lines*78]));

    if(adjointFlag)
    {
      char loadFileAdjoint[MAXLINE];
      write_adjoint_load_file(iter, loadFileAdjoint);
      fprintf(matFile, "*include\n%s\n", loadFileAdjoint);
    }

    fprintf(matFile, "*END");

    fclose(matFile);

    std::cout << "  Mat file " << matFileName << " successfully written\n";

  }

  void compliance_opt::write_adjoint_load_file(const unsigned int iter, char *loadFileName)
  {
    char iter_char[8];
    sprintf(iter_char, "%u", iter);
    strcpy(iter_dir, ".");

    strcpy(loadFileName, "ls_adjoint_load_");
    strcat(loadFileName, iter_char);

    strcat(loadFileName, ".k");

    FILE* loadFile;
    loadFile = fopen(loadFileName, "w");
    if(loadFile == NULL)
    {
      std::cout << "Error opening file, Exiting.\n" << std::endl;
      exit(1);
    }

    fprintf(loadFile, "*keyword\n");
    fprintf(loadFile, "*DEFINE_CURVE_TITLE\ntheForce\n");
    fprintf(loadFile, "$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint\n");
    fprintf(loadFile, "         1         0       1.0       1.0       0.0       0.0         0         0\n");
    for(unsigned int i = 0; i < timeHistory->numSteps; i ++)
    {
      unsigned int nextIndex = timeHistory->numSteps - i - 1;
      double nextVal = -dispPower*pow(timeHistory->u[nextIndex][forcedDofs[0]], dispPower-1.0);
      double nextTime = T_max - timeHistory->t[nextIndex];
      fprintf(loadFile, "    %16f    %16f\n", nextTime, nextVal);
    }

    fprintf(loadFile, "*END");
    fclose(loadFile);
  }


  void compliance_opt::read_displacements()
  {
    std::fill(u.begin(), u.end(), 0.0);


    char dispFileName[MAXLINE];
    strcpy(dispFileName, iter_dir);
    strcat(dispFileName, "/nodout");

    FILE* dispFile;
    dispFile = fopen(dispFileName, "r");
    if(dispFile == NULL)
    {
      std::cout << "Error opening displacement file, Exiting.\n" << std::endl;
      exit(1);
    }

    int dummy;
    std::vector<int> tmp(DIM);
    char nextLine[MAXLINE];
    unsigned int count = 0;
    getNextDataLine(dispFile, nextLine, MAXLINE, &dummy, false);
    while(dummy == 0)
    {
      if(strncmp(nextLine, " nodal point", 12) == 0)
      {
        count ++;
        if(count == 3)
          break;
      }
      getNextDataLine(dispFile, nextLine, MAXLINE, &dummy, false);
    }
    if(dummy == 1)
    {
      std::cout << "Error reading displacement file. Exiting. \n" << std::endl;
      exit(-1);
    }


    unsigned int valuesWritten = 0;
    char tmpch[13];
    for(unsigned int i = 0; i < N_nodes; i ++)
    {
      getNextDataLine(dispFile, nextLine, MAXLINE, &dummy, false);
      if(dummy == 1)
      {
        std::cout << "Error reading displacement file. Exiting. \n" << std::endl;
        exit(-1);
      }

      unsigned int nodeNum;
      valuesWritten = sscanf(nextLine, "%u", &nodeNum);
      if(valuesWritten != 1)
      {
        std::cout << "Error reading displacement file. Exiting. \n" << std::endl;
        exit(-1);
      }


      for(unsigned int k = 0; k < DIM; k ++)
      {
        if(node_ind_dofs[nodeNum - 1][k] == -1)
          continue;

        strncpy(tmpch, &(nextLine[10 + k*12]), 12);
        tmpch[12] = '\0';
        valuesWritten = sscanf(tmpch, "%f", &(u[node_ind_dofs[nodeNum - 1][k]]));
        if(valuesWritten != 1)
        {
          std::cout << "Error reading displacement file. Exiting. \n" << std::endl;
          exit(-1);
        }
      }
    }
    fclose(dispFile);
    std::cout << "  Sucessfully read in displacement file. \n";

  }

  void compliance_opt::set_element_stiffness()
  {
    char stiffFileName[MAXLINE];
    if(!dynamicFlag)
      strcpy(stiffFileName, iter_dir);
    else
      strcpy(stiffFileName, dynamic_init_dir);

    strcat(stiffFileName, "/ElmStfMtx_0001_000");

    FILE* stiffFile;
    stiffFile = fopen(stiffFileName, "r");
    if(stiffFile == NULL)
    {
      std::cout << "Error opening stiffness file, Exiting.\n" << std::endl;
      exit(1);
    }

    unsigned int num_nodes = pow(2, DIM);
    char nextLine[MAXLINE];
    int fileEndFlag;
    unsigned int numStiffnessLines = (DIM*num_nodes*(DIM*num_nodes + 1))/8;
    for(unsigned int i = 0; i < N; i ++)
    {
      std::vector<unsigned int> nextLocalInds(num_nodes);
      std::vector< float > nextStiffnessMat(numStiffnessLines*4);
      unsigned int elmID;

      getNextDataLine(stiffFile, nextLine, MAXLINE, &fileEndFlag, false);
      while(fileEndFlag == false)
      {
        if(strncmp(nextLine, "Element Id", 10) == 0)
          break;

        getNextDataLine(stiffFile, nextLine, MAXLINE, &fileEndFlag, false);
      }
      if(fileEndFlag == true)
      {
        std::cout << "Error 1 reading stiffness file, Exiting.\n" << std::endl;
        exit(1);
      }

      char*  elmID_str = &nextLine[11];
      unsigned int valuesWritten = sscanf(elmID_str, "%u", &elmID);
      if(valuesWritten != 1)
      {
        std::cout << "Error 1 reading stiffness file, Exiting.\n" << std::endl;
        exit(1);
      }


      for(unsigned int k = 0; k < 4; k ++)
        getNextDataLine(stiffFile, nextLine, MAXLINE, &fileEndFlag, false);

      char* nodeNums = &nextLine[12];

      char* tokPtr;
      for(unsigned int k = 0; k < num_nodes; k ++)
      {
        if(k == 0)
          tokPtr = strtok(nodeNums, " ");
        else
          tokPtr = strtok(NULL, " ");

        valuesWritten = sscanf(tokPtr, "%u", &(nextLocalInds[k]));
        if(valuesWritten != 1)
        {
          std::cout << "Error 2 reading stiffness file. Exiting \n" << std::endl;
          exit(-1);
        }
      }

      // line that just says "Element Matrix"
      getNextDataLine(stiffFile, nextLine, MAXLINE, &fileEndFlag, true);
      if(fileEndFlag == true)
      {
        std::cout << "Error 3 reading stiffness file, Exiting.\n" << std::endl;
        exit(1);
      }

      for(unsigned int k = 0; k < numStiffnessLines; k ++)
      {
        getNextDataLine(stiffFile, nextLine, MAXLINE, &fileEndFlag, true);
        if(fileEndFlag == true)
        {
          std::cout << "Error 4 reading stiffness file, Exiting.\n" << std::endl;
          exit(1);
        }

        unsigned int valuesWritten = sscanf(nextLine, "%f %f %f %f",
                        &nextStiffnessMat[k*4],     &nextStiffnessMat[k*4 + 1],
                        &nextStiffnessMat[k*4 + 2], &nextStiffnessMat[k*4 + 3]);
        if(valuesWritten != 4)
        {
          std::cout << "Error 5 reading stiffness file. Exiting \n" << std::endl;
          exit(-1);
        }

      }

      bool foundFlag = false;
      for(unsigned int k = 0; k < N; k ++)
      {
        if(ElmDatas[k].elmID == elmID)
        {
         foundFlag = true;
         ElmDatas[k].set_elm_stiffness(nextLocalInds, node_ind_dofs,
            nextStiffnessMat);
        }
      }
      if(foundFlag == false)
      {
        std::cout << "Error 6 reading stiffness file. Exiting \n" << std::endl;
        exit(-1);
      }

    }

  }


  void compliance_opt::read_mass_matrix()
    {
      char massFileName[MAXLINE];
      strcpy(massFileName, iter_dir);
      strcat(massFileName, "/ElmMssMtx_0001_000");

      FILE* massFile;
      massFile = fopen(massFileName, "r");
      if(massFile == NULL)
      {
        std::cout << "Error opening mass file, Exiting.\n" << std::endl;
        exit(1);
      }

      char nextLine[MAXLINE];
      int fileEndFlag;
      getNextDataLine(massFile, nextLine, MAXLINE, &fileEndFlag, false);

      while(fileEndFlag == false)
      {
        if(strncmp(nextLine, "Lumped Nodal", 12) == 0)
          break;

        getNextDataLine(massFile, nextLine, MAXLINE, &fileEndFlag, false);
      }
      if(fileEndFlag == true)
      {
        std::cout << "Error 1 reading mass file, Exiting.\n" << std::endl;
        exit(1);
      }

      int nextNodeNum = 0;
      double nextMass;
      std::vector<Eigen::Triplet<double> > coefficients;

      for(unsigned int i = 0; i < N_nodes; i ++)
      {
        getNextDataLine(massFile, nextLine, MAXLINE, &fileEndFlag, false);

        unsigned int valuesWritten = sscanf(nextLine, "%i %lf", &nextNodeNum, &nextMass);
        if(valuesWritten != 2)
        {
          std::cout << "Error 2 reading mass file, Exiting.\n" << std::endl;
          exit(1);
        }

        for(unsigned int k = 0; k < DIM; k ++)
        {
          if(node_ind_dofs[nextNodeNum - 1][k] == -1)
            continue;
          else
            coefficients.push_back(
                Eigen::Triplet<double> (node_ind_dofs[nextNodeNum - 1][k],
                                        node_ind_dofs[nextNodeNum - 1][k],
                                        nextMass));
        }

      }

      fclose(massFile);

      M.resize(N_dofs, N_dofs);
      M.setFromTriplets(coefficients.begin(), coefficients.end());
      std::cout << "  Sucessfully set the mass matrix. \n";

    }


  float compliance_opt::compute_objective()
  {

    compliance = 0.0;
    for (unsigned int i = 0; i < N; i++)
    {
      compliance += ElmDatas[i].compute_elm_compliance(u);
    }

    return compliance;
  }

  float compliance_opt::compute_dynamic_objective()
  {
    compliance = 0.0;
    std::vector<double> compliance_history(timeHistory->numSteps, 0.0);
    for(unsigned int k = 0; k < timeHistory->numSteps; k ++)
    {
//      for (unsigned int i = 0; i < K.outerSize(); ++i)
//        for (Eigen::SparseMatrix<double>::InnerIterator it(K,i); it; ++it)
//        {
//          compliance_history[k] += (timeHistory->u[k][it.row()])*(it.value())*(timeHistory->u[k][it.col()]);
//
//        }
//
//      for (unsigned int i = 0; i < M.outerSize(); ++i)
//        for (Eigen::SparseMatrix<double>::InnerIterator it(M,i); it; ++it)
//        {
//          compliance_history[k] += (timeHistory->a[k][it.col()])*(it.value())*(timeHistory->u[k][it.col()]);
//
//        }
      for(unsigned int i = 0; i < forcedDofs.size(); i ++)
        compliance_history[k] += timeHistory->u[k][forcedDofs[i]]*loadingFunction( -timeHistory->t[k], F_max, T_max);

//      std::cout << loadingFunction( -timeHistory->t[k], F_max, T_max) <<  "   " << timeHistory->u[k][forcedDofs[0]]  << " " << forcedDofs.size() <<std::endl;
//
//      std::cout << "     " << compliance_history[k] << std::endl;
    }

    for(unsigned int i = 0; i < timeHistory->numSteps - 1; i ++)
    {
      compliance += 0.5*(compliance_history[i] + compliance_history[i + 1])
                          *(timeHistory->t[i + 1] - timeHistory->t[i]);
    }

//    std::cout << "we in dis bitch \n" << timeHistory->numSteps << std::endl;;


    // check dynamics...
//    double dynamicsWorks;
//    for(unsigned int k = 0; k < timeHistory->numSteps; k ++)
//    {
//      double dynamicsWorks_k = 0.0;
//      std::vector<double> nextStuff(N_dofs, 0.0);
//
//      for (unsigned int i = 0; i < M.outerSize(); ++i)
//        for (Eigen::SparseMatrix<double>::InnerIterator it(M,i); it; ++it)
//        {
//          nextStuff[it.col()] += (timeHistory->a[k][it.col()])*(it.value());
//        }
//
//      for (unsigned int i = 0; i < K.outerSize(); ++i)
//        for (Eigen::SparseMatrix<double>::InnerIterator it(K,i); it; ++it)
//        {
//          nextStuff[it.row()] += it.value()*timeHistory->u[k][it.col()];
//        }
//
//      nextStuff[forcedDofs[0]] -= loadingFunction( -timeHistory->t[k], F_max, T_max);
////      std::cout << " time : " << timeHistory->t[k] << "   load : " << loadingFunction( -timeHistory->t[k], F_max, T_max) <<
////          "     " << timeHistory->a[k][forcedDofs[0]] << std::endl;
////      std::cout << "     : " << nextStuff[forcedDofs[0]] << "   " << nextStuff[forcedDofs[0] + 1] << std::endl;
//
//      for(unsigned int l = 0; l <N_dofs; l ++)
//        dynamicsWorks_k += fabs(nextStuff[l]);
//
////      std::cout << k <<  "  " << dynamicsWorks_k << std::endl;
//      dynamicsWorks += dynamicsWorks_k;
//    }
////
//    std::cout << "dynamics dont work : " << dynamicsWorks << std::endl;


    for (unsigned int i = 0; i < M.outerSize(); ++i)
      for (Eigen::SparseMatrix<double>::InnerIterator it(M,i); it; ++it)
      {
        if(it.col() != it.row())
          std::cout << "   MASS MATRIX PROBLEM : " << it.row() << "  " << it.col();

      }


    return compliance;
  }

  void compliance_opt::compute_sensitivities()
  {

    std::fill(sensitivities.begin(), sensitivities.end(), 0.0);
    for(unsigned int i = 0; i < N; i ++)
      sensitivities[i] = ElmDatas[i].compute_sensitivity(u);

  }

  void compliance_opt::set_stiffness_matrix()
  {
    K.setZero();
    K.data().squeeze();

    std::vector<Eigen::Triplet<double> > coefficients;
    for(unsigned int i =0 ; i < N; i ++)
    {
      double nextScaleFactor = 1.0*ElmDatas[i].E_;
      for(unsigned int a = 0; a < ElmDatas[i].num_dofs; a++)
      {
        if(ElmDatas[i].localToGlobal[a] == -1)
          continue;

        for(unsigned int b = 0; b < ElmDatas[i].num_dofs; b++)
        {
          if(ElmDatas[i].localToGlobal[b] == -1)
            continue;

          coefficients.push_back(
              Eigen::Triplet<double> (ElmDatas[i].localToGlobal[a],
                                      ElmDatas[i].localToGlobal[b],
                                      nextScaleFactor*ElmDatas[i].elmStiffnessMat[a][b]));
        }

      }
    }
    K.setFromTriplets(coefficients.begin(), coefficients.end());

    std::cout << "  Set global stiffness matrix for iteration " << internal_iter << std::endl;
  }

  void compliance_opt::set_mass_matrix()
  {
    M.setZero();
    M.data().squeeze();
    std::vector<Eigen::Triplet<double> > coefficients;

    for(unsigned int i =0 ; i < N; i ++)
    {
      double nextMass = 1.0*(density*ElmDatas[i].volume_*ElmDatas[i].rho_/(pow(2.0, DIM)));
      for(unsigned int a = 0; a < ElmDatas[i].num_dofs; a++)
      {
        if(ElmDatas[i].localToGlobal[a] == -1)
          continue;

        coefficients.push_back(
            Eigen::Triplet<double> (ElmDatas[i].localToGlobal[a],
                                    ElmDatas[i].localToGlobal[a],
                                    nextMass));
      }
    }

    M.setFromTriplets(coefficients.begin(), coefficients.end());
    std::cout << "  Set global mass matrix for iteration " << internal_iter << std::endl;

  }


  void compliance_opt::compute_hessian()
  {
    // first lets invert the stiffness matrix
    Eigen::UmfPackLU<Eigen::SparseMatrix<double> > lu_of_A;
    lu_of_A.compute(K);
    if(!(lu_of_A.info() == Eigen::Success))
    {
      // decomposiiton failed
      std::cout << "Matrix inversion failed. Exiting. \n";
      exit(-1);
    }

    for(unsigned int i = 0; i < N; i ++)
    {
      hessian.push_back(std::vector<double> (N, 0.0));
      Eigen::VectorXd rhs_i(N_dofs);
      Eigen::VectorXd u_tilde_i(N_dofs);

      for(unsigned int k = 0; k < N_dofs; k ++)
        rhs_i[k] = 0.0;

      for(unsigned int a = 0; a < ElmDatas[i].num_dofs; a++)
        for(unsigned int b = 0; b < ElmDatas[i].num_dofs; b ++)
        {
          if(ElmDatas[i].localToGlobal[a] < 0 || ElmDatas[i].localToGlobal[b] < 0)
            continue;
          else
            rhs_i[ElmDatas[i].localToGlobal[a]] +=
                ElmDatas[i].elmStiffnessMat[a][b]*u[ElmDatas[i].localToGlobal[b]];
        }

      rhs_i *= 2.0*(ElmDatas[i].penal_/ElmDatas[i].rho_);

      // solve adjoint problem
      u_tilde_i = lu_of_A.solve(rhs_i);
      if(!(lu_of_A.info() == Eigen::Success))
      {
        // decomposiiton failed
        exit(-1);
      }

      for(unsigned int j = 0; j < N; j ++)
      {
        double multFactor = (ElmDatas[j].penal_/ElmDatas[j].rho_);


        for(unsigned int a = 0; a < ElmDatas[i].num_dofs; a++)
          for(unsigned int b = 0; b < ElmDatas[i].num_dofs; b ++)
          {
            if(ElmDatas[j].localToGlobal[a] < 0 || ElmDatas[j].localToGlobal[b] < 0)
              continue;
            else
              hessian[i][j] += multFactor*
                  u_tilde_i[ElmDatas[j].localToGlobal[a]]*
                    ElmDatas[j].elmStiffnessMat[a][b]*u[ElmDatas[j].localToGlobal[b]];
          }

        if(i == j)
        {
          double multFactor_2 = multFactor*((ElmDatas[j].penal_- 1.0)/ElmDatas[j].rho_);
          for(unsigned int a = 0; a < ElmDatas[j].num_dofs; a++)
            for(unsigned int b = 0; b < ElmDatas[j].num_dofs; b ++)
            {
              if(ElmDatas[j].localToGlobal[a] < 0 || ElmDatas[j].localToGlobal[b] < 0)
                continue;
              else
                hessian[i][j] -= multFactor_2*
                    u[ElmDatas[j].localToGlobal[a]]*
                      ElmDatas[j].elmStiffnessMat[a][b]*u[ElmDatas[i].localToGlobal[b]];
            }
        }
      }
    }
  }

  void compliance_opt::filter_sensitivities()
  {
    std::fill(sensitivities_filtered.begin(), sensitivities_filtered.end(), 0.0);

    for(unsigned int i = 0; i < N ; i ++)
    {
      float H_tot = R_filter;
      sensitivities_filtered[i] += R_filter*ElmDatas[i].rho_*ElmDatas[i].sensitivity_;
      for(unsigned int k = 0; k < ((ElmDatas[i]).neighborList).size(); k++)
      {
        float nextH = R_filter - ((ElmDatas[i]).neighborList)[k].second;
        H_tot += nextH;
        sensitivities_filtered[i] += nextH*(ElmDatas[i].neighborList)[k].first->rho_*(ElmDatas[i].neighborList)[k].first->sensitivity_;
      }
      sensitivities_filtered[i] *= (1.0/(ElmDatas[i].rho_*H_tot));
//      if(sensitivities_filtered[i] > 0.0)
//        sensitivities_filtered[i] = 0.0;
    }

    if(!MMA_flag)
    {
      std::vector<double>::iterator maxIt;
      maxIt = std::max_element(sensitivities_filtered.begin(), sensitivities_filtered.end());
      double maxVal =  *maxIt;
      if(maxVal > 0.0)
        for(unsigned int beta = 0; beta < N; beta++)
          sensitivities_filtered[beta] -= maxVal;
    }

  }

  void compliance_opt::integrate_dynamic_lambda()
  {

    clock_t timer = clock();

    timeHistory->clear_lambda_data();

    std::vector< void *> params;
    params.push_back((void *) &K);
    params.push_back((void *) &M);
    params.push_back((void *) &F_max);
    params.push_back((void *) &T_max);
    params.push_back((void *) &forcedDofs);

    gsl_odeiv2_system sys = {func, jac, 2*N_dofs, &params};

    gsl_odeiv2_driver * d =
      gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                    1e-6, 1e-6, 0.0);
    int i;
    double t = -T_max, t1 = 0.0;
    double *lambda;
    lambda = (double *) calloc(2*N_dofs, sizeof(double));

    timeHistory_adjoint->lambda.resize(numIntegrationSteps + 1);
    timeHistory_adjoint->lambda[numIntegrationSteps] = (double *) calloc(2*N_dofs, sizeof(double));
    timeHistory_adjoint->lambda_time.resize(numIntegrationSteps + 1, 0.0);
    for (i = 1; i <= numIntegrationSteps; i++)
      {
        double ti = -T_max + i * T_max / (1.0*numIntegrationSteps);
        int status = gsl_odeiv2_driver_apply (d, &t, ti, lambda);

        if (status != GSL_SUCCESS)
          {
            printf ("error, return value=%d\n", status);
            break;
          }

        timeHistory_adjoint->lambda_time[numIntegrationSteps - i] = -ti;
        timeHistory_adjoint->lambda[numIntegrationSteps - i] = (double *) calloc(2*N_dofs, sizeof(double));
        for(unsigned int j = 0; j < 2*N_dofs; j ++)
          timeHistory_adjoint->lambda[numIntegrationSteps - i][j] = lambda[j];

      }

    free(lambda);
    gsl_odeiv2_driver_free (d);

    timer = clock() - timer;

    std::cout << "  Sucessfully Integrated the lambda ODE. Duration :  " << ((float)timer)/CLOCKS_PER_SEC << " seconds.\n";

  }

  void compliance_opt::compute_dynamic_dhdrho()
  {

    dh_drho.clear();
    dh_drho_indicies.clear();
    dh_drho_indicies.resize(N, std::vector<unsigned int> (0));
    for(unsigned int timeIndex = 0; timeIndex < timeHistory->numSteps; timeIndex ++)
    {
      double dimFactor = pow(2.0, DIM);
      double dimFactor_inv = 1.0/dimFactor;

      for( unsigned int beta = 0; beta < N; beta++)
      {
        unsigned int elmID = beta; // lmDatas[beta].elmID - 1;
        double next_dK_factor = ElmDatas[beta].dE_drho_; // penal*pow(ElmDatas[beta].rho_, penal - 1.0);
        for(unsigned int i = 0; i < ElmDatas[beta].localToGlobal.size(); i ++)
        {
          int nextGlobalIndex = ElmDatas[beta].localToGlobal[i];
          if( nextGlobalIndex == -1)
            continue;

          std::pair<unsigned int, unsigned int> nextPair = std::make_pair(nextGlobalIndex, elmID);

          if(timeIndex == 0)
          {
            dh_drho.insert(std::make_pair(nextPair,
                                          std::vector<double> (timeHistory->numSteps, 0.0)));
            dh_drho_indicies[elmID].push_back(nextGlobalIndex);
          }
          dh_drho[nextPair][timeIndex] +=  density*dimFactor_inv*ElmDatas[beta].volume_
                                            *timeHistory->a[timeIndex][nextGlobalIndex];

          for(unsigned int j = 0; j < ElmDatas[beta].localToGlobal.size(); j ++)
          {
            int nextGlobalIndex_2 = ElmDatas[beta].localToGlobal[j];

            if( nextGlobalIndex_2 == -1)
              continue;

            dh_drho[nextPair][timeIndex] += next_dK_factor*
                ElmDatas[beta].elmStiffnessMat[i][j]*timeHistory->u[timeIndex][nextGlobalIndex_2];

          }
        }
      }
    }
  }


  void compliance_opt::integrate_dynamic_sensitivity()
  {
    // zero out sensitivities
    std::fill(sensitivities.begin(), sensitivities.end(), 0.0);

    unsigned int numLambdaSteps = timeHistory_adjoint->lambda_time.size();
    for(unsigned int beta = 0; beta < N; beta++)
    {
      std::vector<double> lambdah_beta( numLambdaSteps, 0.0);

      for(unsigned int i = 0; i < dh_drho_indicies[beta].size(); i ++)
      {
        std::pair<unsigned int, unsigned int> nextPair = std::make_pair(dh_drho_indicies[beta][i], beta);

        for(unsigned int j = 0; j < numLambdaSteps; j ++)
        {
          lambdah_beta[j] += timeHistory_adjoint->lambda[j][dh_drho_indicies[beta][i] + N_dofs]*
                interpolate(timeHistory->t, dh_drho[nextPair], timeHistory_adjoint->lambda_time[j], true);
        }
      }
      for(unsigned int j = 0; j < numLambdaSteps-1; j ++)
      {
        sensitivities[beta] += 0.5*(lambdah_beta[j] + lambdah_beta[j+1])*
            (timeHistory_adjoint->lambda_time[j+1] - timeHistory_adjoint->lambda_time[j]);
      }

      ElmDatas[beta].sensitivity_ = sensitivities[beta];

    }
  }

  void compliance_opt::update_rho()
  {
    double l1 = 0.0;
    double l2 = 100000.0;
    double move = 0.2;
    std::vector<float> rho_new(N, 0.0);
    float V = 0.0;
    while ((l2-l1) > 1e-7)
    {
      double lmid = 0.5*(l2+l1);
      for(unsigned int i = 0; i < N; i ++)
        rho_new[i] = std::max((double) rhoMin, std::max(rho[i] - move,
             std::min( 1.0,std::min( (rho[i] + move),
                 (rho[i]*sqrt(-sensitivities_filtered[i]/lmid))))));

      V = 0.0;
      for(unsigned int i = 0; i < N; i ++)
        V += rho_new[i]*ElmDatas[i].volume_;

      if ((V - volFrac*V_tot) > 0.0)
        l1 = lmid;
      else
        l2 = lmid;
    }
    std::cout << "  Volume : " << V << std::endl;

    change = 0.0;
    for(unsigned int i = 0; i < N; i ++)
    {
      if(fabs(rho[i] - rho_new[i]) > change)
        change = fabs(rho[i] - rho_new[i]);
    }
    rho = rho_new;
  }

  void compliance_opt::update_rho_MMA()
  {
    double *xmin = new double[N];
    double *xmax = new double[N];

    double f = compliance;
    double *df = new double[N];

    double g = 0.0;
    double *dg = new double[N];

    double *x = new double[N];

    for(unsigned int i = 0; i < N; i++)
    {
      df[i] = sensitivities_filtered[i];

      xmin[i] = ((rho[i] - 0.2) > rhoMin ? rho[i] - 0.2 : rhoMin);
      xmax[i] = ((rho[i] + 0.2) < 1.0 ? rho[i] + 0.2 : 1.0);

      g += rho[i]*ElmDatas[i].volume_;
      dg[i] = ElmDatas[i].volume_;

      x[i] = rho[i];
    }

    g -= (1.000001)*volFrac*V_tot;

    MMA->Update(x,df,&g,dg,xmin,xmax);

    change = 0.0;
    double V = 0.0;
    for(unsigned int i = 0; i < N; i ++)
    {
      if(fabs(rho[i] - x[i]) > change)
        change = fabs(rho[i] - x[i]);

      rho[i] = x[i];

      V += x[i]*ElmDatas[i].volume_;
    }

    delete [] xmin;
    delete [] xmax;
    delete [] df;
    delete [] dg;
    delete [] x;

    std::cout << "  Volume : " << V << std::endl;

  }

  void compliance_opt::update_element_rhos()
  {
    for(unsigned int i = 0; i < ElmDatas.size(); i ++)
      ElmDatas[i].update_rho(rho[i], k1, k2);
  }

  void compliance_opt::update_element_neighbors()
  {
    for(unsigned int i = 0; i < N; i ++)
    {
      for(unsigned int k = i+1; k < N; k ++)
        if(k == i)
          continue;
        else
        {
          float next_dist = ElmDatas[i].dist(&(ElmDatas[k]));
          if(next_dist < R_filter)
          {
            ElmDatas[i].add_neighbor(&(ElmDatas[k]), next_dist);
            ElmDatas[k].add_neighbor(&(ElmDatas[i]), next_dist);
          }
        }
    }
  }

  void compliance_opt::update_element_vols()
  {
    V_tot = 0.0;

    char d3hspFileName[MAXLINE];
    strcpy(d3hspFileName, iter_dir);
    strcat(d3hspFileName, "/d3hsp");

    FILE* d3hspFile;
    d3hspFile = fopen(d3hspFileName, "r");
    if(d3hspFile == NULL)
    {
      std::cout << "Error opening d3hsp file, Exiting.\n" << std::endl;
      exit(1);
    }

    int fileEndFlag;
    std::vector<int> tmp(DIM);
    char nextLine[MAXLINE];
    getNextDataLine(d3hspFile, nextLine, MAXLINE, &fileEndFlag, false);
    while(fileEndFlag == 0)
    {
      if(strncmp(nextLine, " summary of mass", 16) == 0)
        break;

      getNextDataLine(d3hspFile, nextLine, MAXLINE, &fileEndFlag, false);
    }
    if(fileEndFlag == 1)
    {
      std::cout << "Error 1 reading d3hsp file. Exiting. \n" << std::endl;
      exit(-1);
    }

    char* tokPtr;
    for(unsigned int i = 0; i < N; i ++)
    {
      getNextDataLine(d3hspFile, nextLine, MAXLINE, &fileEndFlag, false);
      if(fileEndFlag == 1)
      {
        std::cout << "Error 2 reading d3hsp file. Exiting. \n" << std::endl;
        exit(-1);
      }
      tokPtr = strtok(nextLine, "=");
      tokPtr = strtok(NULL, "=");
      unsigned int eid;
      float nextVol;
      unsigned int valuesWritten = sscanf(tokPtr, "%u", &eid);
      if(valuesWritten != 1)
      {
        std::cout << "Error 3 reading d3hsp file. Exiting. \n" << std::endl;
        exit(-1);
      }
      tokPtr = strtok(NULL, "=");
      valuesWritten = sscanf(tokPtr, "%f", &nextVol);
      if(valuesWritten != 1)
      {
        std::cout << "Error 4 reading d3hsp file. Exiting. \n" << std::endl;
        exit(-1);
      }

      bool foundFlag = false;
      for(unsigned int k = i; k < N+i; k ++)
      {
        if(ElmDatas[k%N].elmID == eid)
        {
          foundFlag = true;
          ElmDatas[k%N].set_vol(nextVol);
          V_tot += nextVol;
        }
      }
      if(foundFlag == false)
      {
        std::cout << "Error 5 reading d3hsp file. Exiting. \n" << std::endl;
        exit(-1);
      }

    }

  }

  void compliance_opt::write_vtk(unsigned int iter)
  {
    char vtkFileName[MAXLINE];
    char iterChar[MAXLINE];
    sprintf(iterChar, "%u", iter);
//    strcpy(vtkFileName, run_dir);
    strcpy(vtkFileName, "../output/");

    struct stat st;
    if (stat(vtkFileName, &st) == -1)
      mkdir(vtkFileName, 0700);
    strcat(vtkFileName, "/iter_");
    strcat(vtkFileName, iterChar);
    strcat(vtkFileName, ".vtk");

    FILE* vtkFile;
    vtkFile = fopen(vtkFileName, "w");
    if(vtkFile == NULL)
    {
      std::cout << "Error opening vtk file : ";
      std::cout << vtkFileName;
      std::cout << ", Exiting.\n" << std::endl;
      exit(1);
    }
    unsigned int numVert = pow(2, DIM);

    fprintf(vtkFile, "# vtk DataFile Version 3.0\n#This file was generated by the ANDY library\nASCII\nDATASET UNSTRUCTURED_GRID\n\n");
    fprintf(vtkFile, "POINTS %u float\n", N*numVert);
    for(unsigned int i = 0; i < N; i ++)
    {
      for(unsigned int j = 0; j < numVert; j++)
      {
        for(unsigned int k = 0; k < DIM; k ++)
          fprintf(vtkFile, "%f ", node_pos[ElmDatas[i].node_indicies[j] - 1][k]);

        if(DIM < 3)
          fprintf(vtkFile, "0");

        fprintf(vtkFile, "\n");
      }
    }
    fprintf(vtkFile, "\nCELLS %lu %lu\n", ElmDatas.size(), (ElmDatas.size() + N*numVert));
    unsigned int count = 0;
    for(unsigned int i = 0; i < N; i ++)
    {
      fprintf(vtkFile, "%u", numVert);
      for(unsigned int k = 0; k < numVert; k ++)
      {
        fprintf(vtkFile, "\t%u", count);
        count ++;
//        fprintf(vtkFile, "\t%u", (ElmDatas[i].node_indicies[k] - 1));
      }
      fprintf(vtkFile, "\n");
    }
    fprintf(vtkFile, "\n");
    unsigned int type = (DIM == 2 ? 9 : 12);
    fprintf(vtkFile, "CELL_TYPES %u\n", N);
    for(unsigned int i = 0; i < N; i++)
      fprintf(vtkFile, " %u", type);
    fprintf(vtkFile, "\n");

    fprintf(vtkFile, "CELL_DATA %u\n", N);
    fprintf(vtkFile, "SCALARS rho float\n");
    fprintf(vtkFile, "LOOKUP_TABLE default\n");
    for(unsigned int i = 0; i < N; i++)
      fprintf(vtkFile, " %f", ElmDatas[i].rho_);
    fprintf(vtkFile, "\n");
    fprintf(vtkFile, "SCALARS sensitivity float\n");
    fprintf(vtkFile, "LOOKUP_TABLE default\n");
    for(unsigned int i = 0; i < N; i++)
      fprintf(vtkFile, " %f", sensitivities[i]);
    fprintf(vtkFile, "\n");

    fclose(vtkFile);
  }


  void compliance_opt::write_dakota_input(char* static_dir)
  {

    struct stat st;
    if (stat(run_dir, &st) == -1)
      mkdir(run_dir, 0700);


    char relative_static_dir[MAXLINE];
    strcpy(relative_static_dir, "../");
    strcat(relative_static_dir, static_dir);

    chdir(run_dir);
    initialize(relative_static_dir);
    iterate(0);
    chdir("..");

    FILE* dakota_input_file = fopen("compliance_op_dakota.in", "w");
    if(dakota_input_file == NULL)
    {
      std::cout << "Error writing dakota file, exiting.\n";
      exit(-1);
    }

    fprintf(dakota_input_file, "method\n");
    fprintf(dakota_input_file, "  optpp_q_newton\n\n");


    fprintf(dakota_input_file, "variables\n  continuous_design  %u\n", N);
    fprintf(dakota_input_file, "  descriptors  ");
    for(unsigned int i = 0; i < N; i ++)
      fprintf(dakota_input_file, "\n    'rho_%u'", i);

    fprintf(dakota_input_file, "\n  initial_point");
    for(unsigned int i = 0; i < N; i ++)
      fprintf(dakota_input_file, "\n    %f", volFrac);

    fprintf(dakota_input_file, "\n  upper_bounds");
    for(unsigned int i = 0; i < N; i ++)
      fprintf(dakota_input_file, "\n    1.0");

    fprintf(dakota_input_file, "\n  lower_bounds");
    for(unsigned int i = 0; i < N; i ++)
      fprintf(dakota_input_file, "\n    %f", rhoMin);

    fprintf(dakota_input_file, "\n  linear_inequality_constraint_matrix");
    for(unsigned int i = 0; i < N; i ++)
      fprintf(dakota_input_file, "\n    %f", ElmDatas[i].volume_);

    fprintf(dakota_input_file, "\n  linear_inequality_upper_bounds\n    %f\n\n", (volFrac*V_tot));

    fprintf(dakota_input_file, "responses\n  descriptors 'compliance'\n  objective_functions 1\n");
    fprintf(dakota_input_file, "  analytic_gradients\n  analytic_hessians\n\n");

    fprintf(dakota_input_file, "interface\n");
    fprintf(dakota_input_file, "  analysis_driver = \"../../build/dakota_runner\"\n");
    fprintf(dakota_input_file, "  fork\n");
    fprintf(dakota_input_file, "  parameters_file = 'params.in'\n");
    fprintf(dakota_input_file, "  results_file = 'results.out'\n");
    fprintf(dakota_input_file, "  file_save\n");
    fprintf(dakota_input_file, "  deactivate active_set_vector\n");
    fprintf(dakota_input_file, "  work_directory\n");
    fprintf(dakota_input_file, "  named 'workdir'\n");
    fprintf(dakota_input_file, "  directory_save directory_tag\n");

    fclose(dakota_input_file);
  }

  unsigned int compliance_opt::read_dakota_param_input()
  {
    FILE* inFile = fopen("params.in", "r");
    if(inFile == NULL)
    {
      std::cout << "Error opening dakota input file, Exiting.\n";
      exit(-1);
    }

    char nextLine[MAXLINE];
    int dummy;
    getNextDataLine(inFile, nextLine, MAXLINE, &dummy, true);

    for (unsigned int i = 0; i < N; i++)
    {
      getNextDataLine(inFile, nextLine, MAXLINE, &dummy, true);
      sscanf(nextLine, "%f", &rho[i]);
    }

    while(strncmp( &(nextLine[strlen(nextLine) - 8]), "eval_id", 7) != 0)
    {
      getNextDataLine(inFile, nextLine, MAXLINE, &dummy, true);
      if(dummy == 1)
      {
        std::cout << "Error reading dakota input file, Exiting.\n";
        exit(-1);
      }
    }
    unsigned int iter;
    sscanf(nextLine, "%u", &iter);

    fclose(inFile);

    return iter;
  }

  void compliance_opt::write_dakota_param_output()
  {
    FILE* outFile = fopen("results.out", "w");
    if(outFile == NULL)
    {
      std::cout << "Error opening dakota output file, Exiting.\n";
      exit(-1);
    }
    fprintf(outFile, "%f compliance\n[", compliance);
    for (unsigned int i = 0; i < N; i++)
    {
      if(i != 0)
        fprintf(outFile, "\n ");

      fprintf(outFile, " %f", sensitivities[i]);
    }
    fprintf(outFile, " ]\n");

    fprintf(outFile, "[[");
    for (unsigned int i = 0; i < N; i++)
    {
      for (unsigned int j = 0; j < N; j++)
      {
        fprintf(outFile, "\n ");
        fprintf(outFile, " %f", hessian[i][j]);

      }
    }
    fprintf(outFile, " ]]\n");


    fclose(outFile);
  }


  void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                          int const maxSize, int* const endOfFileFlag, bool trimFlag)
  {
    *endOfFileFlag = 0;
    do
    {
      if(fgets(nextLinePtr, maxSize, filePtr) == NULL)
      {
        *endOfFileFlag = 1;
        break;
      }
      if(trimFlag)
      {
        while ((nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t') ||
               (nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r' ))
        {
          nextLinePtr = (nextLinePtr + 1);
        }
      }
    }
    while ((strncmp("~", nextLinePtr, 1) == 0)
           || (strncmp("#", nextLinePtr, 1) == 0)
           || (strncmp("$", nextLinePtr, 1) == 0)
           || (strlen(nextLinePtr) == 0));
  }

}


#endif
