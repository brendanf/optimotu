// SPDX-FileCopyrightText: 2025 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_SEQUENCEVIEW_R_H_INCLUDED
#define OPTIMOTU_SEQUENCEVIEW_R_H_INCLUDED

#ifdef OPTIMOTU_R

#include "SequenceView.h"
#include <Rcpp.h>

SequenceSet sequence_views_from_r(const Rcpp::CharacterVector &seq);

#endif // OPTIMOTU_R

#endif // OPTIMOTU_SEQUENCEVIEW_R_H_INCLUDED
