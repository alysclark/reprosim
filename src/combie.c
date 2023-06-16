//libraries//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

#include "placentastdlib.h"

// merged Main + plac_flow, basically. We couple the detailed placenta model with the fetal circulation model.
// I will only comment the bits of the code that are not directly copied from the above two. For the rest, see the
// 'decoupled' codes first to understand what's going on.

int run(void)
{

    int num_Anodes = 63608;
    int num_Aelts = 63607;
    int num_term = 31803;

    // placenta vars now have '_pl' to distinguish them from fetal
    float C_pl[2*num_Anodes];
    double V_pl[2*num_Anodes];
    float Rvec[2*num_Aelts + num_term];
    int templinks[2*(2*num_Aelts + num_term)];
    int tempnode_elts[3*(2*num_Anodes)];
    int i = 0;

    double Qvec[2*num_Aelts + num_term];
    double P_pl[2*num_Anodes];
    float time_spent=0.0;
    int k;
    // these two are the flows that couple the placenta to the adjacent fetal compartments
    double Q_in;
    double Q_out;
    double flow;
    char filename[100];
    char ending[] = ".txt";
    float R_join = 0.285; // resistance given to links that tie the placenta to the adjacent fetal compartments.
    // I took the resistance of Pennati on either side of the placenta, subtracted the total R of the placenta
    // (computed in 'placenta.m') and divided by 2.
    float U_v = 40;
    float U_a = 3;
    float E_sys = 3;
    float E_dia = 0.3;
    float R_v = 0.08;

    float T_c = 0.43;
    float DeltaT = T_c/4;
    float T_v_s = T_c/2;
    // float T_v_d = T_c - T_v_s;
    float T_a_s = T_c/4;
    // float T_a_d = T_c - T_a_s;

    float dt = 1E-7; // tiny step needed for accurate evolution
    float rough_dt = 1E-2; // writing out more densely kills hard drive space too quickly
    int rough_steps = 101; // number of rough steps
    int num_tsteps = (rough_steps-1)*round(rough_dt/dt)+1; // number of fine steps
    double t;
    int ind;
    // note: previously time-dependent functions are now evaluated one t-step at a time to save memory space (not necessary to have them as a vector)
    double A_a;
    double A_v;
    double Ua;
    double Uv;
    double E;
    double td;
    float C[23]={0};
    float L[23][23]={0};
    float R[23][23]={0};
    float K[23][23]={0};
    float D[23][23]={0};

    // Note: absolutely crucial to have main time-dep vars as 'double' precision, otherwise calc goes off
    double P[23] = {21.8, 21.9, 2.1, 3.0, 44.0, 44.4, 43.1, 41.4, 40.8, 40.3, 43.1, 11.0, 42.1, 32.4, 4.7, 19.3, 5.1, 11.1, 32.1, 4.3, 22.4, 6.7, 10.2};
    double V[23]={0};
    double P_new[23]={0};
    double Q[23][23]={0};
    double Q_new[23][23]={0};
    int indy;
    int jones;
    double DP;
    double disc;
    double Q1;
    double blah;
    int step;

    FILE *file;

    //directly from matlab (computed therein)//
    int indy1[44] = {14, 19, 11, 7, 15, 6, 8, 7, 9, 16, 17, 18, 8, 20, 22, 11, 3, 10, 13, 12, 14, 2, 13, 15, 6, 14, 8, 17, 19, 21, 8, 16, 8, 19, 2, 16, 18, 22, 9, 21, 16, 20, 9, 19};
    int jones1[44] = {2, 2, 3, 6, 6, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 10, 11, 11, 12, 13, 13, 14, 14, 14, 15, 15, 16, 16, 16, 16, 17, 17, 18, 18, 19, 19, 19, 19, 20, 20, 21, 21, 22, 22};
    int indy2[6] = {4, 5, 19, 0, 1, 3};
    int jones2[6] = {0, 1, 3, 4, 5, 19};
    int indy3[2] = {21, 19};
    int jones3[2] = {19, 21};
    int indy4[12] = {2, 3, 0, 1, 10, 6, 5, 12, 10, 4, 7, 6};
    int jones4[12] = {0, 1, 2, 3, 4, 5, 6, 6, 7, 10, 10, 12};
    int indy5[8] = {2, 3, 14, 19, 11, 19, 0, 1};
    int jones5[8] = {0, 1, 2, 2, 3, 3, 4, 5};

    file = fopen("C.txt", "r");
    for (i=0; i<2*num_Anodes; i++)
        fscanf(file, "%e", &C_pl[i]);
    fclose(file);

    file = fopen("V.txt", "r");
    for (i=0; i<2*num_Anodes; i++)
        fscanf(file, "%le", &V_pl[i]);
    fclose(file);

    file = fopen("Rvec.txt", "r");
    for (i=0; i<2*num_Aelts + num_term; i++)
        fscanf(file, "%e", &Rvec[i]);
    fclose(file);

    file = fopen("templinks.txt", "r");
    for (i=0; i<2*(2*num_Aelts + num_term); i++)
        fscanf(file, "%d", &templinks[i]);
    fclose(file);

    file = fopen("tempnode_elts.txt", "r");
    for (i=0; i<3*(2*num_Anodes); i++)
        fscanf(file, "%d", &tempnode_elts[i]);
    fclose(file);

    for (i = 0; i<(2*num_Anodes); i++)
        P_pl[i] = V_pl[i]/C_pl[i];

    for (i = 0; i<(2*num_Aelts+num_term); i++)
        Qvec[i] = (P_pl[templinks[i]-1] - P_pl[templinks[2*num_Aelts + num_term+i]-1])/Rvec[i];

    step = 0;
    t = dt*step;

    td = fmod(t,T_c);
    if (td<0) {td = td + T_c;}
    if (td<=T_a_s) {A_a = sin(pi/T_a_s*td);}
    else {A_a = 0.0;}

    td = fmod((t-DeltaT),T_c);
    if (td<0) {td = td + T_c;}
    if (td<=T_v_s) {A_v = sin(pi/T_v_s*td);}
    else {A_v = 0.0;}

    Ua = U_a*A_a;
    Uv = U_v*A_v;

    E = E_dia + E_sys*A_v;

    C[0] = 1;
    C[1] = 1;
    C[2] = 2;
    C[3] = 1;
    C[4] = 0.08;
    C[5] = 0.05;
    C[6] = 0.08;
    C[7] = 0.07;
    C[8] = 0.04;
    C[9] = 0.05;
    C[10] = 0.08;
    C[11] = 0.4;
    C[12] = 0.01;
    C[13] = 0.3;
    C[14] = 1;
    C[15] = 0.85;
    C[16] = 3;
    C[17] = 0.25;
    C[18] = 0.02;
    C[19] = 0.6;
    C[20] = 1.5;
    C[21] = 0.3;
    C[22] = 4;

    L[1][3] = 0.0016; L[3][1] = 0.0016;
    L[0][2] = 0.0016; L[2][0] = 0.0016;
    L[7][10] = 0.006; L[10][7] = 0.006;
    L[6][12] = 0.08; L[12][6] = 0.08;
    L[5][6] = 0.002; L[6][5] = 0.002;
    L[4][10] = 0.002; L[10][4] = 0.002;

    R[5][6] = 0.12; R[6][5] = 0.12;
    R[6][7] = 0.4; R[7][6] = 0.4;
    R[7][8] = 0.04; R[8][7] = 0.04;
    R[8][9] = 0.06; R[9][8] = 0.06;
    R[4][10] = 0.07; R[10][4] = 0.07;
    R[10][11] = 13.5; R[11][10] = 13.5;
    R[7][10] = 0.01; R[10][7] = 0.01;
    R[6][12] = 0.3; R[12][6] = 0.3;
    R[12][13] = 3; R[13][12] = 3;
    R[13][14] = 8.5; R[14][13] = 8.5;
    R[6][15] = 8; R[15][6] = 8;
    R[14][15] = 4.9; R[15][14] = 4.9;
    R[8][16] = 81; R[16][8] = 81;
    R[8][17] = 34; R[17][8] = 34;
    R[16][17] = 7; R[17][16] = 7;
    R[8][18] = 3.5; R[18][8] = 3.5;
    R[18][19] = 14; R[19][18] = 14;
    R[9][20] = 3.9; R[20][9] = 3.9;
    R[20][21] = 3.4; R[21][20] = 3.4;
    R[9][22] = 3.5; R[22][9] = 3.5;
    R[19][22] = 0.6; R[22][19] = 0.6;
    R[21][16] = 0.5; R[16][21] = 0.5;
    R[16][19] = 0.16; R[19][16] = 0.16;
    R[19][21] = 1.3; R[21][19] = 1.3;
    R[2][14] = 0.2; R[14][2] = 0.2;
    R[2][19] = 0.12; R[19][2] = 0.12;
    R[3][11] = 2; R[11][3] = 2;

    K[1][3] = 0.002; K[3][1] = 0.002;
    K[0][2] = 0.002; K[2][0] = 0.002;
    K[0][4] = 0.001; K[4][0] = 0.001;
    K[1][5] = 0.001; K[5][1] = 0.001;
    K[3][19] = 0.4; K[19][3] = 0.4;
    K[7][10] = 0.009; K[10][7] = 0.009;
    K[19][21] = 0.26; K[21][19] = 0.26;

    D[0][4] = 1;
    D[1][5] = 1;
    D[2][0] = 1;
    D[3][1] = 1;
    D[11][3] = 1;
    D[19][3] = 1;
    D[19][2] = 1;
    D[14][2] = 1;

    // IC for P is given with variable initialisation as it's more efficient //
    for (ind = 0; ind<23; ind++)
        V[ind] = P[ind]*C[ind];

    blah = 0;
    for (i = 0; i<2*num_Anodes; i++)
        blah = blah + V_pl[i];
    V[20] = blah; // track total placenta volume

    P[0] = Uv + E*V[0];
    P[1] = Uv + E*V[1];
    P[2] = Ua + V[2]/C[2];
    P[3] = Ua + V[3]/C[3];

    for (indy = 0; indy<23; indy++)
    {
        for (jones = 0; jones<23; jones++)
        {

            if (R[indy][jones]!=0 && K[indy][jones]==0) {
                DP = P[indy] - P[jones];
                Q[indy][jones] = DP/R[indy][jones];
            }

            if (R[indy][jones]==0 && K[indy][jones]!=0) {
                DP = P[indy]-P[jones];
                Q[indy][jones] = sqrt(fabs(DP)/K[indy][jones])*copysign(1.0,DP);
                if ((indy==19 && jones==3) || (indy==3 && jones==19)) {
                    Q[indy][jones] = pow(fabs(DP)/K[indy][jones],1.6)*copysign(1.0,DP);
                }
            }

            if (R[indy][jones]!=0 && K[indy][jones]!=0) {
                DP = P[indy] - P[jones];
                disc = pow(R[indy][jones],2) + 4*K[indy][jones]*fabs(DP);
                if (disc<0) {
                    Q[indy][jones] = 0;
                }
                else {
                    Q1 = (-R[indy][jones]+sqrt(disc))/(2*K[indy][jones]);
                    if (Q1>0) {
                        Q[indy][jones] = Q1*copysign(1.0,DP);
                    }
                    else {
                        Q[indy][jones] = 0;
                    }
                }
            }

        }
    }



    for (ind = 0; ind<8; ind++) {
        if (Q[indy5[ind]][jones5[ind]]<0) {
            Q[indy5[ind]][jones5[ind]] = 0;
            Q[jones5[ind]][indy5[ind]] = 0;
        }
    }

    // coupling flows btw fetal compartments and the placenta
    Q_in = (P[9]-P_pl[0])/R_join;
    Q_out = (P_pl[num_Anodes]-P[21])/R_join;

    // put these into Q by hand to make sure the rest of the fetal model evolves correctly
    Q[9][20] = Q_in; Q[20][9] = -Q_in;
    Q[21][20] = -Q_out; Q[20][21] = Q_out;


    for (ind = 0; ind<23; ind++)
    {P_new[ind] = P[ind];}


    for (indy = 0; indy<23; indy++) {
        for (jones = 0; jones<23; jones++) {
            Q_new[indy][jones] = Q[indy][jones];
        }
    }

    // write out initial point
    chdir("data");
    itoa(1,filename,10);
    strcat(filename,ending);

    file=fopen(filename,"w");
    for(i=0;i<2*num_Anodes;i++)
        fprintf(file,"%le\n",V_pl[i]);
    fclose(file);

    chdir("V");
    file=fopen(filename,"w");
    for(i=0;i<23;i++)
        fprintf(file,"%le\n",V[i]);
    fclose(file);
    chdir("..");

    chdir("P");
    file=fopen(filename,"w");
    for(i=0;i<23;i++)
        fprintf(file,"%le\n",P[i]);
    fclose(file);
    chdir("..");

    chdir("Q");
    file=fopen(filename,"w");
    for(i=0;i<23;i++) {
        for(k=0;k<23;k++) {
            fprintf(file,"%le\n",Q[i][k]);
        }
    }
    fclose(file);
    chdir("..");

    chdir("..");

    clock_t begin = clock();

    // loop over time
    for (step = 1; step<num_tsteps; step++) {

        t = dt*step;
        // evaluate t-dep functions
        td = fmod(t,T_c);
        if (td<0) {td = td + T_c;}
        if (td<=T_a_s) {A_a = sin(pi/T_a_s*td);}
        else {A_a = 0.0;}

        td = fmod((t-DeltaT),T_c);
        if (td<0) {td = td + T_c;}
        if (td<=T_v_s) {A_v = sin(pi/T_v_s*td);}
        else {A_v = 0.0;}

        Ua = U_a*A_a;
        Uv = U_v*A_v;

        E = E_dia + E_sys*A_v;
        // evolve placenta first
        for (jones = 0; jones<2*num_Anodes; jones++)
        {
            flow = 0;
            for (k = 0; k<3; k++)
            {
                if ((templinks[tempnode_elts[jones+k*2*num_Anodes]-1]==jones+1) && (tempnode_elts[jones+k*2*num_Anodes]!=0)) {
                    flow = flow - Qvec[tempnode_elts[jones+k*2*num_Anodes]-1];}
                else if ((templinks[tempnode_elts[jones+k*2*num_Anodes]-1+2*num_Aelts + num_term]==jones+1) && (tempnode_elts[jones+k*2*num_Anodes]!=0)) {
                    flow = flow + Qvec[tempnode_elts[jones+k*2*num_Anodes]-1];}
            }
            // note: noise filter used in plac_flow is not onlt unhelpful here, it is harmful!
            if (jones==0) {flow = flow + Q_in;}
            else if (jones==num_Anodes) {flow = flow - Q_out;}
            // observe in/out flows (coupling)
            V_pl[jones] = V_pl[jones] + dt*flow;
        }

        for (i = 0; i<(2*num_Anodes); i++)
            P_pl[i] = V_pl[i]/C_pl[i];

        for (i = 0; i<(2*num_Aelts+num_term); i++)
            Qvec[i] = (P_pl[templinks[i]-1] - P_pl[templinks[i+2*num_Aelts + num_term]-1])/Rvec[i];
        // now evolve fetal model
        blah = 0;
        for (indy = 0; indy<23; indy++)
            blah = blah + Q[indy][0];

        P_new[0] = Uv + E*V[0] + R_v*blah;

        blah = 0;
        for (indy = 0; indy<23; indy++)
            blah = blah + Q[indy][1];

        P_new[1] = Uv + E*V[1] + R_v*blah;
        P_new[2] = Ua + V[2]/C[2];
        P_new[3] = Ua + V[3]/C[3];


        for (jones = 4; jones<23; jones++)
        {
            blah = 0;
            for (indy = 0; indy<23; indy++)
                blah = blah + Q[indy][jones];

            P_new[jones] = P[jones] + dt/C[jones]*blah;
        }

        for (ind = 0; ind<44; ind++) {
            DP = P[indy1[ind]] - P[jones1[ind]];
            Q_new[indy1[ind]][jones1[ind]] = DP/R[indy1[ind]][jones1[ind]];
        }

        for (ind = 0; ind<6; ind++) {
            DP = P[indy2[ind]] - P[jones2[ind]];
            Q_new[indy2[ind]][jones2[ind]] = sqrt(fabs(DP)/K[indy2[ind]][jones2[ind]])*copysign(1.0,DP);
            if ((indy2[ind]==19 && jones2[ind]==3) || (indy2[ind]==3 && jones2[ind]==19)) {
                Q_new[indy2[ind]][jones2[ind]] = pow(fabs(DP)/K[indy2[ind]][jones2[ind]],1.6)*copysign(1.0,DP);
            }
        }

        for (ind = 0; ind<2; ind++) {
            DP = P[indy3[ind]]-P[jones3[ind]];
            disc = pow(R[indy3[ind]][jones3[ind]],2) + 4*K[indy3[ind]][jones3[ind]]*fabs(DP);
            if (disc<0) {
                Q_new[indy3[ind]][jones3[ind]] = 0;}
            else {
                Q1 = (-R[indy3[ind]][jones3[ind]]+sqrt(disc))/(2*K[indy3[ind]][jones3[ind]]);
                if (Q1>0) {
                    Q_new[indy3[ind]][jones3[ind]] = Q1*copysign(1.0,DP);}
                else {
                    Q_new[indy3[ind]][jones3[ind]] = 0;
                }
            }
        }


        for (ind = 0; ind<12; ind++) {
            DP = P[indy4[ind]]-P[jones4[ind]];
            Q_new[indy4[ind]][jones4[ind]] = Q[indy4[ind]][jones4[ind]] + (DP-R[indy4[ind]][jones4[ind]]*Q[indy4[ind]][jones4[ind]]-K[indy4[ind]][jones4[ind]]*pow(Q[indy4[ind]][jones4[ind]],2)*copysign(1.0,Q[indy4[ind]][jones4[ind]]))/L[indy4[ind]][jones4[ind]]*dt;
        }

        for (ind = 0; ind<8; ind++) {
            if (Q_new[indy5[ind]][jones5[ind]]<0) {
                Q_new[indy5[ind]][jones5[ind]] = 0;
                Q_new[jones5[ind]][indy5[ind]] = 0;
            }
        }
        // coupling through the flows
        Q_in = (P_new[9]-P_pl[0])/R_join;
        Q_out = (P_pl[num_Anodes]-P_new[21])/R_join;

        Q_new[9][20] = Q_in; Q_new[20][9] = -Q_in;
        Q_new[21][20] = -Q_out; Q_new[20][21] = Q_out;


        for (jones = 0; jones<23; jones++) {
            blah = 0;
            for (indy = 0; indy<23; indy++)
                blah = blah + Q[indy][jones];

            V[jones] = V[jones] + dt*blah;
        }

        blah = 0;
        for (i = 0; i<2*num_Anodes; i++)
            blah = blah + V_pl[i];
        V[20] = blah;

        for (jones = 0; jones<23; jones++) {
            P[jones] = P_new[jones];
            for (indy = 0; indy<23; indy++)
                Q[jones][indy] = Q_new[jones][indy];
        }

        // write out data
        if (step % 100000 == 0){ //the 100000 is hard-coded "round(rough_dt/dt)"
            chdir("data");
            itoa(round(t/rough_dt)+1,filename,10);
            strcat(filename,ending);

            file=fopen(filename,"w");
            for(i=0;i<2*num_Anodes;i++)
                fprintf(file,"%le\n",V_pl[i]);
            fclose(file);

            chdir("V");
            file=fopen(filename,"w");
            for(i=0;i<23;i++)
                fprintf(file,"%le\n",V[i]);
            fclose(file);
            chdir("..");

            chdir("P");
            file=fopen(filename,"w");
            for(i=0;i<23;i++)
                fprintf(file,"%le\n",P[i]);
            fclose(file);
            chdir("..");

            chdir("Q");
            file=fopen(filename,"w");
            for(i=0;i<23;i++) {
                for(k=0;k<23;k++) {
                    fprintf(file,"%le\n",Q[i][k]);
                }
            }
            fclose(file);
            chdir("..");

            chdir("..");
        }

    }

    clock_t end = clock();

    time_spent += (float)(end - begin) / CLOCKS_PER_SEC;

    printf("The elapsed time is %f seconds", time_spent);

    return 0;
}

int main(void)
{
    return run();
}
