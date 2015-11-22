#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <GetOpt.h>

#include "codonString.h"

const char* helpStr = 
    "This program simulates an evolving viral population where the epidemiology\n"
    "is described by SIRS dynamics and the molecular evolution is given by a codon model.\n"
    "The program outputs a fastA file containing viral DNA sequences sampled\n"
    "from the population\n\n"
    "Usage: epiPop -o <outputfile> [options]\n\n"
    "Options:\n"
    "  -h \t\t\t print help message\n"
    "  -g <numeric > 0>\t rate of host immune decay per timestep [0.1]\n"
    "  -c <0<numeric<1>\t strength of immune protection [0.1]\n"
    "  -b <numeric > 0>\t rate of viral infection [0.5]\n"
    "  -d <numeric > 0>\t rate of infection recovery [0.1]\n"
    "  -N <integer > 0>\t initial viral population size [10]\n"
    "  -H <integer > 0>\t host population size [10000]\n"
    "  -L <integer > 0>\t length of viral codon sequence [num antigenic sites]\n"
    "  -l <integer>\t\t number of antigenic codons in viral sequence [2]\n"
    "  -u <numeric>\t\t mutation rate per codon per timestep [1e-5]\n"
    "  -k <numeric>\t\t transition-transversion rate of mutational model [3.0]\n"
    "  -a <numeric>\t\t cost of antigenic mutation [1e-2] \n"
    "  -s <numeric>\t\t cost of non-antigenic mutation [1e-3]\n"
    "  -T <numeric > 0>\t max running time of simulation [10000]\n"
    "  -S <numeric > 0>\t interval between sampling [1000]\n"
    "  -Z <integer>\t\t number of sequences sampled [10]\n"
    "  -t <numeric > 0>\t size of time-step [1.0]\n\n";

void seed_time (void);

int main (int argc, char **argv)
{
    // Default simulation parameters
    char filename[100] = "";
    double Tmax = 10000;                        // length of simulation
    double Tsample = 1000;                      // time of first sample (and intervals)
    int sampleSize = 10;                        // number of seqs sampled
    double dt = 1.0;                            // time-step


    // Default epidemiological parameters
    int hostPopSize = 10000;                     
    int N0 = 10;                                // initial viral pop size
    double gamma = 0.1;                         // decay of immunity
    double beta = 0.5;                          // transmission rate
    double delta = 0.1;                         // recovery rate

    // Default sequence parameters
    double kappa = 3.0;
    double mut_rate = 1e-5;                     
    int Ld = 0;                                 
    int Lb = 2;                              
    double sb = 0.01;
    double sd = 1e-3;
    double sigma = 1.0;

    // User-specified arguments
    int opt_char;
    while ((opt_char = getopt(argc, argv, "ho:T:S:Z:t:H:N:g:b:d:L:u:k:l:a:s:c:")) != -1)
    {
        switch (opt_char)
        {
            case 'h':
                printf("%s\n", helpStr);
                exit(0);
                break;
            case 'o':
                strcpy(filename, optarg);
                break;
            case 'T':
                if ((Tmax = strtod(optarg, NULL)) <= 0.0)
                {
                    fprintf (stderr, "Invalid -T parameter: %s\n", optarg);
                    exit(1);
                }
                break;

            case 'S':
                if ((Tsample = strtod(optarg, NULL)) <= 0.0)
                {
                    fprintf (stderr, "Invalid -S parameter: %s\n", optarg);
                    exit(1);
                }
                break;
            case 'Z':
                sampleSize = atoi(optarg);
                break;
            case 't':
                if ((dt = strtod(optarg, NULL)) <= 0.0)
                {
                    fprintf (stderr, "Invalid -t parameter: %s\n", optarg);
                    exit(1);
                }
                break; 
           case 'N':
                if ((N0 = atoi(optarg)) <= 0)
                {
                    fprintf (stderr, "Invalid -N parameter: %s\n", optarg);
                    exit(1);
                }
                break;
           case 'H':
                if ((hostPopSize = atoi(optarg)) <= 0)
                {
                    fprintf (stderr, "Invalid -H parameter: %s\n", optarg);
                    exit(1);
                }
                break;

            case 'g':
                if ((gamma = strtod(optarg, NULL)) < 0.0)
                {
                    fprintf (stderr, "Invalid -g parameter: %s\n", optarg);
                    exit(1);
                }
                break;
            case 'b':
                if ((beta = strtod(optarg, NULL)) <= 0.0)
                {
                    fprintf (stderr, "Invalid -b parameter: %s\n", optarg);
                    exit(1);
                }
                break;
            case 'd':
                if ((delta = strtod(optarg, NULL)) <= 0.0)
                {
                    fprintf (stderr, "Invalid -d parameter: %s\n", optarg);
                    exit(1);
                }
                break;
            case 'L':
                if ((Ld = atoi(optarg)) < 0)
                {
                    fprintf (stderr, "Invalid -L parameter: %s\n", optarg);
                    exit(1);
                }
                break;
            case 'u':
                if ((mut_rate = strtod(optarg, NULL)) < 0.0)
                {
                    fprintf (stderr, "Invalid -u parameter: %s\n", optarg);
                    exit(1);
                }
                break;
            case 'k':
                if ((kappa = strtod(optarg, NULL)) < 0.0)
                {
                    fprintf (stderr, "Invalid -k parameter: %s\n", optarg);
                    exit(1);
                }
                break;
            case 'l':
                if ((Lb = atoi(optarg)) < 0)
                {
                    fprintf (stderr, "Invalid -l parameter: %s\n", optarg);
                    exit(1);
                }
                break;
            case 'a':
                sb = strtod(optarg, NULL);
                break;
            case 's':
                sd = strtod(optarg, NULL);
                break;
            case 'c':
                if ((sigma = strtod(optarg, NULL)) < 0.0)
                {
                    fprintf (stderr, "Invalid -c parameter: %s\n", optarg);
                    exit(1);
                }
                break;
            case '?':
                fprintf (stderr, "Unrecognized argument\n");
                exit(1);
        }
    }

    if (filename[0] == '\0')
    {
        fprintf (stderr, "Output filename (option -o) must be specified\n");
        exit(1);
    }
    FILE *outfile;
    if ((outfile = fopen(filename, "w")) == NULL)
    {
        fprintf (stderr, "Cannot open %s\n", filename);
        exit(1);
    }

    seed_time();
    codonString::init_trans_matrix (kappa);
    codonString viralPop(N0, Lb, Ld, sb, sd, sigma);

    double t = 0.0;
    int N = 0;
    while ((N = viralPop.popSize()) > 0 && t < Tmax)
    {
        viralPop.mutate(mut_rate * dt);
        viralPop.transmit(beta * dt, hostPopSize);
        viralPop.recover(rbinom(N, delta * dt));

        if (gamma > 0)
            viralPop.immuneDecay(gamma * dt);

        if (fmod(t, Tsample) < 1e-4)
            viralPop.printSample(outfile, t, sampleSize);
        t += dt;
    }
    fclose(outfile);
}
