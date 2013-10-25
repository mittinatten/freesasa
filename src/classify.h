#ifndef SASALIB_CLASSIFY_H
#define SASALIB_CLASSIFY_H

double classify_radius(const char *res_name, const char *atom_name);

int classify_class(const char *res_name, const char *atom_name);

const char* classify_class2str(int class);

int classify_residue(const char *res_name, const char *atom_name);

const char* classify_residue2str(int res);

int classify_element(const char *atom_name);

const char* classify_element2str(int element);

double classify_element_radius(int element);

int classify_oons(const char *res_name, const char *atom_name);

const char* classify_oons2str(int oons_type);

double classify_oons_radius(int oons_type);

int classify_is_aminoacid(const char *res_name);

int classify_is_nucleotide(const char *res_name);

int classify_is_polar(int class);

int classify_is_apolar(int class);

#endif
