#ifndef MY_COMPLIANCE_OPT_C
#define MY_COMPLIANCE_OPT_C

#include "compliance_opt.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include <cmath>
#include <unistd.h>


namespace Comp_Op
{
/******************** ElmData member functions **********************/

  ElmData::ElmData(unsigned int dim, double pen)
  :
  dim_(dim),
  penal_(pen)
  {
    num_vertex = pow(2, dim);
    num_dofs = num_vertex*dim;

    elmStiffnessMat.resize(num_dofs, std::vector<float> (num_dofs));
    localToGlobal.resize(num_dofs);
  }

  void ElmData::set_elm_stiffness(std::vector<unsigned int> &local_inds,
                                  std::vector<std::vector<int> > &node_ind_dofs,
                                  std::vector<float> &stiff_values)
  {
    unsigned int count = 0;
    for(unsigned int i = 0; i < num_vertex; i ++)
      for(unsigned int k = 0; k < dim_; k ++)
      {
        localToGlobal[count] = node_ind_dofs[local_inds[i]][k];
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

  float ElmData::compute_sensitivity(const std::vector<float> &u, float rho_e)
  {
    float sensitivity = 0.0;
    for(unsigned int i = 0; i < num_dofs; i ++)
      for(unsigned int k = 0; k < num_dofs; k ++)
      {
        if(localToGlobal[i] == -1 || localToGlobal[k] == -1)
          continue;
        else
          sensitivity += u[localToGlobal[i]]*elmStiffnessMat[i][k]*u[localToGlobal[k]];
      }

    float multFactor = -penal_/rho_e;
    sensitivity *= multFactor;

    return sensitivity;
  }

  /******************** compliance_opt member functions **********************/

  compliance_opt::compliance_opt()
  {
    volFrac = 0.5;
    rhoMin = 0.0;
    penal = 3.0;

    density = 1.0;
    poisson = 0.3;
    N = 200;
    N_nodes = 231;
    N_dofs = 440;

    firstFlag = true;

    node_ind_dofs.resize(N_nodes, std::vector<int>(2));
    u.resize(N_dofs);
    sensitivities.resize(N);
    ElmDatas.resize(N, ElmData(DIM, penal));

    rho.resize(N);
    std::fill(rho.begin(), rho.end(), volFrac);

  }

  void compliance_opt::write_element_file()
  {
    unsigned int eid, pid, n1, n2, n3, n4, n5, n6, n7, n8;

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

    getNextDataLine(fid, nextLine,MAXLINE, &fileEndFlag);
    valuesWritten = sscanf(nextLine, "%u %u %u %u %u %u %u %u %u %u",
                            &eid, &pid, &n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8);
    if (valuesWritten != 10)
    {
      fileReadErrorFlag = true;
    }

    while(fileEndFlag == false && fileReadErrorFlag == false)
    {
      fprintf(fid_w, "%8u%8u%8u%8u%8u%8u%8u%8u%8u%8u\n",
        eid, eid, n1, n2, n3, n4, n5, n6, n7, n8);

      getNextDataLine(fid, nextLine,MAXLINE, &fileEndFlag);
      valuesWritten = sscanf(nextLine, "%u %u %u %u %u %u %u %u %u %u",
                            &eid, &pid, &n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8);
      if (valuesWritten != 10)
      {
        fileReadErrorFlag = true;
        break;
      }
    }
    fprintf(fid_w, "*END");

    fclose(fid);
    fclose(fid_w);

    if (fileReadErrorFlag)
    {
      std::cout << "Error reading file, Exiting.\n" << std::endl;
      exit(1);
    }
    else
     std::cout << "Successfully wrote element file" << std::endl;

    char tmp[MAXLINE];
    getcwd(ls_static_dir, MAXLINE);

    std::cout << "Element File Sucessfully Written\n";


    chdir("-");

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
//      fprintf(fid_w, "%10u%10u%10u%10u%10u%10u%10u%10u\n\n",
//                       i+1, 1, i+1, 0,0,0,0,0);

      fprintf(fid_w, "%10u%10u%10u\n",
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

    std::cout << "Part File Sucessfully Written\n";

    chdir("-");
  }

  void compliance_opt::run_ls_dyna(const unsigned int iter)
  {
    chdir(run_dir);
    char matFileName[MAXLINE];
    write_mat_file(iter, matFileName);

//    char command[MAXLINE];
//
//    sprintf(command, "cd %s;", iter_dir);
//    strcat(command, "module load bhatta/ls-dyna-double-smp/r1010;");
//
//    char runCommand[MAXLINE];
//    sprintf(runCommand,
//        "/groups/bhatta/software/ls-dyna_smp_d_r1010_x64_redhat5_ifort160/ls-dyna_smp_d_r1010_x64_redhat5_ifort160 i=%s >>ResultOutput.txt",
//        matFileName);
//
//    strcat(command, runCommand);
//
//    system(command);

  }

  void compliance_opt::write_mat_file(const unsigned int iter, char *matFileName)
  {
    char iter_char[8];
    sprintf(iter_char, "%u", iter);
    strcpy(iter_dir, run_dir);
    strcat(iter_dir, "/iter_");
    strcat(iter_dir, iter_char);

    struct stat st;
    if (stat(iter_dir, &st) == -1)
      mkdir(iter_dir, 0700);

    chdir(iter_dir);

    strcpy(matFileName, "ls_mat_");
    strcat(matFileName, iter_char);
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
      fprintf(matFile, "%10u %9g %9g %9g %9g %9g %9u\n", i+1, density, rho[i], poisson,
                  dummy, dummy, dummy_int);

    }
    fprintf(matFile, "*include_path_relative\n../../keyword_static\n");
    fprintf(matFile, "*include\n");
    char ls_info_path[MAXLINE];
    sprintf(ls_info_path, "ls_info.k");
    unsigned int lines = std::floor((1.0*strlen(ls_info_path))/78.0);
    for(unsigned int i = 0; i < lines; i ++)
    {
      char nextLine[78];
      strncpy(nextLine, &(ls_info_path[i*78]), 78);
      fprintf(matFile, "%s +\n", nextLine);
    }
    fprintf(matFile, "%s\n", &(ls_info_path[lines*78]));

    fprintf(matFile, "*END");

    fclose(matFile);

    chdir("-");

  }

  void compliance_opt::read_node_inds()
  {
    char nodeFileName[MAXLINE];
    strcpy(nodeFileName, iter_dir);
    strcat(nodeFileName, "/Node_Data_0001_001");

    FILE* nodeFile;
    nodeFile = fopen(nodeFileName, "r");
    if(nodeFile == NULL)
    {
      std::cout << "Error opening node file, Exiting.\n" << std::endl;
      exit(1);
    }

    int dummy;
    std::vector<int> tmp(DIM);
    char nextLine[MAXLINE];
    for(unsigned int i = 0; i < 4; i ++)
      getNextDataLine(nodeFile, nextLine,MAXLINE, &dummy);

    unsigned int valuesWritten;
    char *tokPtr;
    for(unsigned int i = 0; i < N_nodes; i ++)
    {
      getNextDataLine(nodeFile, nextLine, MAXLINE, &dummy);
      tokPtr = strtok(nextLine, " ");
      for(unsigned int k = 0; k < DIM; k ++)
      {
        tokPtr = strtok(NULL, " ");
        valuesWritten = sscanf(tokPtr, "%d", &(tmp[k]));
        if(valuesWritten != 1)
        {
          std::cout << "Error reading node file. Exiting \n" << std::endl;
          exit(-1);
        }
        node_ind_dofs[i][k] = tmp[k] - 1;
      }
    }

    fclose(nodeFile);
    std::cout << "Sucessfully read node file.\n";
  }

  void compliance_opt::read_displacements()
  {
    if(firstFlag)
    {
      read_node_inds();
      firstFlag = false;
    }

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

      for(unsigned int k = 0; k < DIM; k ++)
      {
        if(node_ind_dofs[i][k] == -1)
          continue;

        strncpy(tmpch, &(nextLine[10 + k*12]), 12);
        tmpch[12] = '\0';
        valuesWritten = sscanf(tmpch, "%f", &(u[node_ind_dofs[i][k]]));
        if(valuesWritten != 1)
        {
          std::cout << "Error reading displacement file. Exiting. \n" << std::endl;
          exit(-1);
        }
      }
    }
    fclose(dispFile);
    std::cout << "Sucessfully read in displacement file. \n";

  }

  void compliance_opt::set_element_stiffness()
  {
    char stiffFileName[MAXLINE];
    strcpy(stiffFileName, iter_dir);
    strcat(stiffFileName, "/ElmStfMtx_0001_000");

    FILE* stiffFile;
    stiffFile = fopen(stiffFileName, "r");
    if(stiffFile == NULL)
    {
      std::cout << "Error opening stiffness file, Exiting.\n" << std::endl;
      exit(1);
    }

    unsigned int num_nodes = pow(2, DIM);
    std::vector<unsigned int> nextLocalInds(num_nodes);
    char nextLine[MAXLINE];
    int fileEndFlag;
    char* nodeNums;
    char* tokPtr;
    unsigned int numStiffnessLines = (DIM*num_nodes*(DIM*num_nodes + 1))/8;
    std::vector< float > nextStiffnessMat(numStiffnessLines*4);
    for(unsigned int i = 0; i < N; i ++)
    {
      getNextDataLine(stiffFile, nextLine, MAXLINE, &fileEndFlag, false);
      while(fileEndFlag == false)
      {
        if(strncmp(nextLine, "  node list", 11) == 0)
          break;

        getNextDataLine(stiffFile, nextLine, MAXLINE, &fileEndFlag, false);
      }
      if(fileEndFlag == true)
      {
        std::cout << "Error 1 reading stiffness file, Exiting.\n" << std::endl;
        exit(1);
      }
      nodeNums = &nextLine[12];

      for(unsigned int k = 0; k < num_nodes; k ++)
      {
        if(k == 0)
          tokPtr = strtok(nodeNums, " ");
        else
          tokPtr = strtok(NULL, " ");

        unsigned int valuesWritten = sscanf(tokPtr, "%u", &(nextLocalInds[k]));
        if(valuesWritten != 1)
        {
          std::cout << "Error 2 reading stiffness file. Exiting \n" << std::endl;
          exit(-1);
        }
        nextLocalInds[k]--;
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

      ElmDatas[i].set_elm_stiffness(nextLocalInds, node_ind_dofs,
          nextStiffnessMat);

    }

  }

  void compliance_opt::compute_sensitivities()
  {
    set_element_stiffness();

    std::fill(sensitivities.begin(), sensitivities.end(), 0.0);
    for(unsigned int i = 0; i < N; i ++)
      sensitivities[i] = ElmDatas[i].compute_sensitivity(u, rho[i]);

  }

  void compliance_opt::getNextDataLine( FILE* const filePtr, char* nextLinePtr,
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
    while ((strncmp("*", nextLinePtr, 1) == 0)
           || (strncmp("#", nextLinePtr, 1) == 0)
           || (strncmp("$", nextLinePtr, 1) == 0)
           || (strlen(nextLinePtr) == 0));
  }

}

#endif
