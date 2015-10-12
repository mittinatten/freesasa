/*
  Copyright Simon Mitternacht 2013-2015.

  This file is part of FreeSASA.

  FreeSASA is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  FreeSASA is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with FreeSASA.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CLASSIFIER_H
#define CLASSIFIER_H

#include "freesasa.h"

/**
    Returns reference to the default classifier, and keeps reference
    count. Every module that calls this function must also call
    freesasa_classifier_default_release() to free resources.
 */
const freesasa_classifier *
freesasa_classifier_default_acquire();

/**
    Frees default classifier (by decreasing reference count), actually
    frees it when count is zero.
 */
void
freesasa_classifier_default_release();

double
freesasa_guess_radius(const char* symbol);

#endif /* CLASSIFIER_H */
