//libraries//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

#include "placentastdlib.h"

// data is written out from combie.c & continue.c as variable vectors at each point in time.
// to plot a variable as a function of time we load them one by one and extract only the var we need.

int run(void)
{
    int num_Anodes = 63608;
    int num_Aelts = 63607;
    int num_term = 31803;

    float C_pl[2*num_Anodes];
    double V_pl[2*num_Anodes];
    // double P_pl[2*num_Anodes];
    float Rvec[2*num_Aelts + num_term];
    // double Qvec[2*num_Aelts + num_term];
    int tempnode_elts[3*(2*num_Anodes)];


    int compt=1; //+num_Anodes; // compt number by MATLAB'S normal standard. ie. starts at 1.
    int compt2=2; //+num_Anodes; // set to zero is you want to examine volume or pressure of a
    // single compt. set to int>0 if you want to look at flow btw the 2.
    int rough_steps = 501; // set to final data point label (in this case, 501)
    int ind;
    int i;
    int k;
    int common;

    // the three vars below will track the compt(s) we want
    double vol[rough_steps];
    double press[rough_steps];
    double vol2[rough_steps];
    double press2[rough_steps];
    double flow[rough_steps];

    char filename[100];
    char ending[] = ".txt";
    int current[3]={0};
    int current2[3]={0};

    FILE *file;

    // load placenta variables for calculations
    file = fopen("C.txt", "r");
    for (i=0; i<2*num_Anodes; i++)
        fscanf(file, "%e", &C_pl[i]);
    fclose(file);

    file = fopen("Rvec.txt", "r");
    for (i=0; i<2*num_Aelts + num_term; i++)
        fscanf(file, "%e", &Rvec[i]);
    fclose(file);

    file = fopen("tempnode_elts.txt", "r");
    for (i=0; i<3*(2*num_Anodes); i++)
        fscanf(file, "%d", &tempnode_elts[i]);
    fclose(file);

    chdir("data");

    // read in volume as a function of time, extract needed compt, and write out final vector

    for (ind = 0; ind<rough_steps; ind++) {
        itoa(ind+1,filename,10);
        strcat(filename,ending);

        file = fopen(filename, "r");
        for (i=0; i<2*num_Anodes; i++) {
            fscanf(file, "%le", &V_pl[i]);
        }
        fclose(file);

        vol[ind] = V_pl[compt-1];
        if (compt2>0) {vol2[ind] = V_pl[compt2-1];}
    }
    chdir("..");

    strcpy(filename,"vol");
    strcat(filename,ending);
    file=fopen(filename,"w");
    for(i=0;i<rough_steps;i++)
        fprintf(file,"%le\n",vol[i]);
    fclose(file);

    // calc pressure and write it out

    for (i = 0; i<rough_steps; i++) {
        press[i] = vol[i]/C_pl[compt-1];
        if (compt2>0) {press2[i] = vol2[i]/C_pl[compt2-1];}
    }

    strcpy(filename,"press");
    strcat(filename,ending);
    file=fopen(filename,"w");
    for(i=0;i<rough_steps;i++)
        fprintf(file,"%le\n",press[i]);
    fclose(file);

    // calc flow and write it out

    if (compt2>0){
        current[0] = tempnode_elts[compt-1];
        current[1] = tempnode_elts[compt-1+2*num_Anodes];
        current[2] = tempnode_elts[compt-1+2*(2*num_Anodes)];
        current2[0] = tempnode_elts[compt2-1];
        current2[1] = tempnode_elts[compt2-1+2*num_Anodes];
        current2[2] = tempnode_elts[compt2-1+2*(2*num_Anodes)];

        for (i=0;i<3;i++){
            for (k=0;k<3;k++){
                if (current[i]==current2[k] && current[i]>0 && current2[k]>0) {
                    common = current[i];
                }
            }
        }

        for (i = 0; i<rough_steps; i++) {
            flow[i] = (press[i] - press2[i])/Rvec[common-1];
        }


        strcpy(filename,"flow");
        strcat(filename,ending);
        file=fopen(filename,"w");
        for(i=0;i<rough_steps;i++)
            fprintf(file,"%le\n",flow[i]);
        fclose(file);
    }


    return 0;
}

int main(void)
{
    return run();
}
