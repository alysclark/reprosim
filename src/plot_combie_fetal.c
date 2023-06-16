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
    int compt=10; // compt number by MATLAB'S normal standard. ie. starts at 1.
    int compt2=23; // set to zero is you want to examine volume or pressure of a
    // single compt. set to int>0 if you want to look at flow btw the 2.
    int rough_steps = 501; // set to final data point label (in this case, 501)
    int ind;
    int i;
    int k;
    // the three vars below will track the compt(s) we want
    double vol[rough_steps];
    double press[rough_steps];
    double flow[rough_steps];

    char filename[100];
    char ending[] = ".txt";

    // to house loaded vars:
    double V[23];
    double P[23];
    double Q[23][23];

    FILE *file;

    // read in volume as a function of time, extract needed compt, and write out final vector

    chdir("data/V");

    for (ind = 0; ind<rough_steps; ind++) {
        itoa(ind+1,filename,10);
        strcat(filename,ending);

        file = fopen(filename, "r");
        for (i=0; i<23; i++)
            fscanf(file, "%le", &V[i]);
        fclose(file);

        vol[ind] = V[compt-1];
    }

    chdir("../..");
    strcpy(filename,"vol");
    strcat(filename,ending);
    file=fopen(filename,"w");
    for(i=0;i<rough_steps;i++)
        fprintf(file,"%le\n",vol[i]);
    fclose(file);

    // read in pressure as a function of time, extract needed compt, and write out final vector

    chdir("data/P");

    for (ind = 0; ind<rough_steps; ind++) {
        itoa(ind+1,filename,10);
        strcat(filename,ending);

        file = fopen(filename, "r");
        for (i=0; i<23; i++)
            fscanf(file, "%le", &P[i]);
        fclose(file);

        press[ind] = P[compt-1];
    }

    chdir("../..");
    strcpy(filename,"press");
    strcat(filename,ending);
    file=fopen(filename,"w");
    for (i=0;i<rough_steps;i++)
        fprintf(file,"%le\n",press[i]);
    fclose(file);

    // read in flow as a function of time, extract needed transfer btw compts, and write out final vector

    if (compt2>0) {
        chdir("data/Q");

        for (ind = 0; ind<rough_steps; ind++) {
            itoa(ind+1,filename,10);
            strcat(filename,ending);

            file = fopen(filename, "r");
            for (i=0;i<23;i++) {
                for (k=0;k<23;k++) {
                    fscanf(file,"%le", &Q[i][k]);
                }
            }
            fclose(file);

            flow[ind] = Q[compt-1][compt2-1];
        }


        chdir("../..");
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
