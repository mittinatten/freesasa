#ifndef TOOLS_H
#define TOOLS_H

// This enum is now only used in tests, and has therefore been removed from freesasa.h
typedef enum {
    FREESASA_APOLAR=0, // Apolar atoms
    FREESASA_POLAR, // Polar atoms
    FREESASA_CLASS_UNKNOWN // neither
} freesasa_class;

int
float_eq(double a, double b, double tolerance);

void
set_fail_after(int freq);

#endif
