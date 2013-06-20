#ifndef CLASSIFIER_H
#define CLASSIFIER_H

typedef struct {
    int nclasses;
    const char **class2str;
    int (*classify)(const char *res_name, const char *atom_name);
} atomclassifier;

#endif
