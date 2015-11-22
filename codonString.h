#ifndef BITSTRING_H
#define BITSTRING_H

#include <vector>
#include "aa_code.h"

class codonString
{
    const int Lb;
    const int Ld;
    const int seq_len;
    const double sb;
    const double sd;
    const double sigmaL;

    double wmean;
    double wsq;

    std::vector<double> RR;
    std::vector<int> ancSeq;

    static double cum_trans[NCODON][NCODON];

    struct seqDiff
    {
        int pos;
        int codon;

        seqDiff(int p, int c) : pos(p), codon(c) {}
        bool operator==(const seqDiff& mut) { return (pos == mut.pos); }
        bool operator==(int p) { return(pos == p); }
    };

    struct seqVar
    {
        int first_nAg;             
        double fitness;
        std::vector<seqDiff> changes;


        void mutateSeq (int pos, int anc_cdn, double s_eff, int Lb);      
        seqVar() : first_nAg(0), fitness(1.0) {}
    };
    std::vector<seqVar> seqList;
 
    int aaj_index (int pos, int aa) const { return (pos*NAA + aa); }
    static int mutateCodon (int old);

  public:

    codonString(int pop_size, int Lb, int Ld, double sb, double sd, double sigma);
    static void init_trans_matrix (double kappa);
    
    void printSample(FILE* pFile, double t, int n) const;
    void immuneDecay(double gamma);
    void mutate(double mu);

    int transmit(double beta, int hostPopSize);
    void recover(int num);
    int popSize() { return (seqList.size()); }

};

#endif /* BITSTRING_H  */
