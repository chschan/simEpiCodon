#ifndef CODON_H_
#define CODON_H_

#define NCODON 64
#define NAA 21
#define STOP '*'

const char codon[NCODON][4] = {"aaa", "aag", "aac", "aat",
                         "aga", "agg", "agc", "agt",
                         "aca", "acg", "acc", "act",
                         "ata", "atg", "atc", "att",
                         "gaa", "gag", "gac", "gat",
                         "gga", "ggg", "ggc", "ggt",
                         "gca", "gcg", "gcc", "gct",
                         "gta", "gtg", "gtc", "gtt",
                         "caa", "cag", "cac", "cat",
                         "cga", "cgg", "cgc", "cgt",
                         "cca", "ccg", "ccc", "cct",
                         "cta", "ctg", "ctc", "ctt",
                         "taa", "tag", "tac", "tat",
                         "tga", "tgg", "tgc", "tgt",
                         "tca", "tcg", "tcc", "tct",
                         "tta", "ttg", "ttc", "ttt"};

const char aa[NAA] = {'K', 'N', 'R', 'S', 'T',
                           'I', 'M', 'E', 'D', 'G',
                           'A', 'V', 'Q', 'H', 'F',
                           'P', 'L', 'Y', 'W', 'C',
                           STOP};

const int cd2aa[NCODON] = {0, 0, 1, 1,
                        2, 2, 3, 3,
                        4, 4, 4, 4,
                        5, 6, 5, 5,
                        7, 7, 8, 8,
                        9, 9, 9, 9,
                        10, 10, 10, 10,
                        11, 11, 11, 11,
                        12, 12, 13, 13,
                        2, 2, 2, 2,
                        15, 15, 15, 15,
                        16, 16, 16, 16,
                        20, 20, 17, 17,
                        20, 18, 19, 19,
                        3, 3, 3, 3,
                        16, 16, 15, 15};

const char AA[NCODON] = {'K', 'K', 'N', 'N',
                   'R', 'R', 'S', 'S',
                   'T', 'T', 'T', 'T',
                   'I', 'M', 'I', 'I',
                   'E', 'E', 'D', 'D',
                   'G', 'G', 'G', 'G',
                   'A', 'A', 'A', 'A',
                   'V', 'V', 'V', 'V',
                   'Q', 'Q', 'H', 'H',
                   'R', 'R', 'R', 'R',
                   'P', 'P', 'P', 'P',
                   'L', 'L', 'L', 'L',
                   STOP, STOP, 'Y', 'Y',
                   STOP, 'W', 'C', 'C',
                   'S', 'S', 'S', 'S',
                   'L', 'L', 'F', 'F'};



#endif /* CODON_H */

