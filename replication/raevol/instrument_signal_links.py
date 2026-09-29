#!/usr/bin/env python3
"""Patch a COPY of the aevol-7 source so that the per-generation statistics of the
best individual also report the regulatory influences exerted by SIGNAL proteins.

R-Aevol adds each signal to the individual's protein list and computes its
affinities like those of any other protein (ExpManager_7.cpp: the condition is
`prot->is_init_ || prot->signal_`), so nb_link and avg_influence include
signal-to-RNA influences. This patch counts them separately and appends two
columns at the very END of every row of stats/stats_simd_best.csv:

  nb_signal_link          number of positive enhancer/operator factors whose
                          protein is a signal (same counting rule as nb_link)
  sum_signal_influence    sum of those factors

All existing columns are computed and written exactly as before, so a run resumed
from a backup must reproduce the original rows column for column; the analysis
checks this before using the new columns. Genome-encoded ("internal") regulation
is then nb_link - nb_signal_link and S - sum_signal_influence.

Usage: instrument_signal_links.py <source root of the COPY>   (e.g. ~/raevol7_instr)
"""
import os, sys

root = sys.argv[1]
h = os.path.join(root, "src/libaevol/7/Stats_7.h")
c = os.path.join(root, "src/libaevol/7/Stats_7.cpp")


def patch(path, edits):
    s = open(path).read()
    for old, new in edits:
        n = s.count(old)
        if n != 1:
            sys.exit("%s: expected 1 match, found %d for:\n%s" % (path, n, old))
        s = s.replace(old, new)
    open(path, "w").write(s)


patch(h, [(
    "    double  av_value_influences_ = 0.0;",
    "    double  av_value_influences_ = 0.0;\n"
    "    int32_t nb_signal_influences_ = 0;     // instrumented: influences of signal proteins\n"
    "    double  sum_signal_influences_ = 0.0;  // instrumented: sum of their factors",
)])

patch(c, [
    ("""  double mean_activator_activity = 0.0;
  double mean_operator_activity = 0.0;
""", """  double mean_activator_activity = 0.0;
  double mean_operator_activity = 0.0;
  int32_t nb_signal_infl = 0;       // instrumented
  double sum_signal_infl = 0.0;     // instrumented
"""),
    ("""          if (affinity.enhancer_factor > 0.0)
          {
            nb_activators++;
            mean_activator_activity += affinity.enhancer_factor;
          }

          if (affinity.operator_factor > 0.0)
          {
            nb_operators++;
            mean_operator_activity += affinity.operator_factor;
          }""", """          bool from_signal = (affinity.protein != nullptr && affinity.protein->signal_);  // instrumented
          if (affinity.enhancer_factor > 0.0)
          {
            nb_activators++;
            mean_activator_activity += affinity.enhancer_factor;
            if (from_signal) { nb_signal_infl++; sum_signal_infl += affinity.enhancer_factor; }
          }

          if (affinity.operator_factor > 0.0)
          {
            nb_operators++;
            mean_operator_activity += affinity.operator_factor;
            if (from_signal) { nb_signal_infl++; sum_signal_infl += affinity.operator_factor; }
          }"""),
    ("""  nb_influences_                 = nb_operating_influences_ + nb_enhancing_influences_;""",
     """  nb_influences_                 = nb_operating_influences_ + nb_enhancing_influences_;
  nb_signal_influences_          = nb_signal_infl;    // instrumented
  sum_signal_influences_         = sum_signal_infl;   // instrumented"""),
    ("""      statfile_best_<< "," << indiv->nb_bases_non_essential_including_nf_genes_;
    }
    statfile_best_<<std::endl;""",
     """      statfile_best_<< "," << indiv->nb_bases_non_essential_including_nf_genes_;
    }
    #ifdef __REGUL
      statfile_best_<< "," << nb_signal_influences_ << "," << sum_signal_influences_;  // instrumented, always last
    #endif
    statfile_best_<<std::endl;"""),
])
print("patched", h, "and", c)
