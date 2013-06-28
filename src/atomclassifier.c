#include <string.h>
#include "atomclassifier.h"

typedef enum {GLY=0, ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO,
	      SER, THR, CYS, CSE, ASN, GLN, HIS, TYR,
	      ASP, GLU, ARG, LYS, PYH, ASX, GLX, residue_unknown} residue_type;

const char *residue_types[] =  
{   "GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO",
    "SER", "THR", "CYS", "CSE", "ASN", "GLN", "HIS", "TYR",
    "ASP", "GLU", "ARG", "LYS", "PYH", "ASX", "GLX", "UNK"};

//typedef enum {residue_hydrophobic,residue_polar} residue_class;

int atomclassifier_residue_classify(const char *res_name, 
				    const char *atom_name);

atomclassifier atomclassifier_residue() 
{
    atomclassifier residue_classifier = 
	{ .nclasses = residue_unknown+1, 
	  .class2str = residue_types, 
	  .classify = &atomclassifier_residue_classify };
    return residue_classifier;
}

int atomclassifier_residue_classify(const char *res_name, 
				    const char *atom_name) 
{
    for (int i = 0; i < residue_unknown+1; ++i) {
	if (! strcmp(res_name,residue_types[i])) return i;
    }
    if (! strcmp(res_name,"SEC")) return CSE;
    if (! strcmp(res_name,"PYL")) return PYH;
    return residue_unknown;
}
