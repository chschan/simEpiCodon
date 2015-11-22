#define MATHLIB_STANDALONE 1
#include <Rmath.h>
#include <stdio.h>
#include <algorithm>
#include <numeric>
#include <stddef.h>

#include "codonString.h"

using namespace std;

double codonString::cum_trans[NCODON][NCODON];

codonString::codonString (int pop_size, int Lb, int Ld, 
    double sb, double sd, double sigma):
    Lb(Lb), Ld(Ld), seq_len(Lb + Ld), sb(sb), sd(sd), sigmaL(sigma/Lb), 
    wmean (1.0), wsq(0.0),
    RR(NAA*Lb, 0.0), ancSeq(seq_len, -1), seqList(pop_size, seqVar())
{
    // Initialise ancestral sequences (assume equal freqs)
    int cdn;
    for (int pos = 0; pos < seq_len; ++pos)
    {
        do
            cdn = runif(0, NCODON);
        while (AA[cdn] == STOP);
        ancSeq[pos] = cdn;
    }
}

void codonString::init_trans_matrix (double kappa)
{
    double tmp = 0;
    for (int i = 0; i < NCODON; ++i)
    {
        for (int j = 0; j < NCODON; ++j)
        {
            if (i == j)
                tmp = 0.0;
            else
            {
                int diffs[3];
                int pos = -1; 
                for (int p = 0; p < 3; ++p)
                {
                    diffs[p] = (codon[i][p] != codon[j][p]);
                    if (diffs[p])
                        pos = p;
                }

                if ((diffs[0] + diffs[1] + diffs[2]) > 1)
                    tmp = 0;
                else
                {
                    tmp = 1;
                    switch (codon[i][pos])
                    {
                        case 'a':
                            if (codon[j][pos] == 'g')
                                tmp = kappa;
                            break;
                        case 'g':
                            if (codon[j][pos] == 'a')
                                tmp = kappa;
                            break;
                        case 'c':
                            if (codon[j][pos] == 't')
                                tmp = kappa;
                            break;
                        case 't':
                            if (codon[j][pos] == 'c')
                                tmp = kappa;
                            break;
                    }
                }
            }

            if (j == 0)
                cum_trans[i][j] = tmp;
            else
                cum_trans[i][j] = cum_trans[i][j-1] + tmp;
        }
    }
}

// Prints ancestral sequence and a random sample (with replacement) of num sequences
void codonString::printSample (FILE* pFile, double t, int num) const
{
    char tmpSeq[(seq_len * 3) + 1];
    tmpSeq[seq_len * 3] = '\0';
    for (int pos = 0; pos < seq_len; ++pos)
        strncpy(&tmpSeq[3*pos], codon[ancSeq[pos]], 3);
    //fprintf(pFile, "> t=0.0, ancSeq\n%s\n", tmpSeq);
    
    for (int i = 0; i < num; ++i)
    {
        vector<seqDiff>::const_iterator it;
        int index = runif(0, seqList.size());
        for (it = seqList[index].changes.begin(); it != seqList[index].changes.end(); ++it)
            strncpy(&tmpSeq[3*(it->pos)], codon[it->codon], 3);
        fprintf (pFile, "> t=%.2f, i=%d, w=%f, I=%zu, wmean=%f, wsq=%f\n%s\n", 
            t, index, seqList[index].fitness, seqList.size(), wmean, wsq, tmpSeq);

        // undo changes
        for (it = seqList[index].changes.begin(); it != seqList[index].changes.end(); ++it)
            strncpy(&tmpSeq[3*(it->pos)], codon[ancSeq[it->pos]], 3);
    }
}


void codonString::immuneDecay(double gamma)
{
    for (int i = 0; i < RR.size(); ++i)
    {
        if (RR[i] > 0)
        {
            RR[i] -= rbinom(RR[i], gamma);
        }
    }
}

void codonString::mutate(double mu)
{
    int num = rpois(mu * seq_len * seqList.size());
    for (int i = 0; i < num; i++)
    {
        int R = (int) runif(0, seqList.size());                 
        int pos = (int) runif(0, seq_len);

        if (pos < Lb)
            seqList[R].mutateSeq(pos, ancSeq[pos], sb, Lb);        
        else
            seqList[R].mutateSeq(pos, ancSeq[pos], sd, Lb);
    }
}

void codonString::seqVar::mutateSeq (int pos, int anc_cdn, double s_eff, int Lb)
{
    int old_cdn = anc_cdn;
    vector<seqDiff>::iterator it = find(changes.begin(), changes.end(), pos);
    if (it != changes.end())
        old_cdn = it->codon;
    int new_cdn = mutateCodon(old_cdn);

    // Update fitness
    if (AA[new_cdn] != AA[anc_cdn])
        fitness = fitness * (1.0 - s_eff);
    if (AA[old_cdn] != AA[anc_cdn])
        fitness = fitness / (1.0 - s_eff);
    if (AA[new_cdn] == STOP)
        fitness = 0.0;
    if (fitness < 0.0)
        fitness = 0.0;

    // Update seqVar
    if (it != changes.end())
    {
        it->codon = new_cdn;

        // check for reversions
        if (it->codon == anc_cdn)
        {
            if (pos < Lb)
            {
                *it = changes[first_nAg-1];
                if (first_nAg < changes.size())
                    changes[first_nAg-1] = changes.back();
                changes.pop_back();
                first_nAg--;
            }
            else
            {
                *it = changes.back();
                changes.pop_back();
                if (first_nAg == changes.size())
                    first_nAg--;
            }
        }
    }
    else
    {
        if (pos < Lb && first_nAg < changes.size())
        {
            changes.push_back(changes[first_nAg]);
            changes[first_nAg] = seqDiff(pos, new_cdn);
            first_nAg++;
        }
        else
        {
            changes.push_back(seqDiff(pos, new_cdn));
            if (pos < Lb || first_nAg == changes.size())
                first_nAg++;
        }
    }

}

int codonString::mutateCodon(int old)
{
    double tmpR = runif(0, cum_trans[old][NCODON-1]);
    int j;
    for (j = 0; j < NCODON; ++j)
    {
        if (old == j)
            continue;
        if (tmpR <= cum_trans[old][j])
            break;
    }
    return (j);
}

// death is a neutral process
void codonString::recover (int num)
{
    for (int i = 0; i < num; ++i)
    {
        if (seqList.size() < 2)
            return;

        // Update memory status
        size_t index = (int) runif(0, seqList.size());
        for (int j = 0; j < Lb; ++j)
            RR[aaj_index(j,cd2aa[ancSeq[j]])]++;

        for (int k = 0; k < seqList[index].first_nAg; ++k)
        {
            int pp = seqList[index].changes[k].pos;
            RR[aaj_index(pp,cd2aa[ancSeq[pp]])]--;
            RR[aaj_index(pp,cd2aa[seqList[index].changes[k].codon])]++;
        }

        if (index != seqList.size() - 1)
            seqList[index] = seqList.back();
         seqList.pop_back();
    }
}


// The fitness of the individual is based on BOTH beneficial and deleterious sites
int codonString::transmit (double beta, int hostPopSize)
{
    if (seqList.size() == 0)
        return (0);
  
    // Get status of population
    int N = seqList.size();
    vector<double> cum_fitness(N, 0);
    cum_fitness[0] = seqList[0].fitness;
    wsq = (seqList[0].fitness*seqList[0].fitness);
    for (int i = 1; i < N; ++i)
    {
        cum_fitness[i] = cum_fitness[i-1] + seqList[i].fitness;
        wsq += (seqList[i].fitness*seqList[i].fitness);
    }
    wmean = cum_fitness[N-1]/(1.0 * N);
    int num = rpois(beta * wmean * N);
    
    // Selectively (based on wD) pick parental sequences
    int real_births = 0;
    vector<double> probs(num, 0.0);
    for (int i = 0; i < num; ++i)
        probs[i] = runif(0, cum_fitness.back());
    sort(probs.begin(), probs.end());

    double ancRR = 0;
    for (int j = 0; j < Lb; ++j)
        ancRR += RR[aaj_index(j, cd2aa[ancSeq[j]])];

    vector<double>::iterator lb = cum_fitness.begin();
    for (int i = 0; i < num; ++i)
    {
        lb = lower_bound(lb, cum_fitness.end(), probs[i]);
        ptrdiff_t index = lb - cum_fitness.begin();

        double tmpRR = ancRR;
        for (int j = 0; j < seqList[index].first_nAg; ++j)
        {
            int pp = seqList[index].changes[j].pos;
            tmpRR -= RR[aaj_index(pp,cd2aa[ancSeq[pp]])];
            tmpRR += RR[aaj_index(pp,cd2aa[seqList[index].changes[j].codon])];
        }
        double pSus = (hostPopSize - N - sigmaL*tmpRR) * 1.0/hostPopSize;

        if (pSus > runif(0,1))
        {
            real_births++;
            seqList.push_back(seqList[index]);
        }
    }
    return (real_births);
}


