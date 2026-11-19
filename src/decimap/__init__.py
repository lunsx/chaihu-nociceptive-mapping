"""DECIMAP: de-circular, evidence-graded computational mapping.

DECIMAP packages the bias-control procedures used in the manuscript. Version 2
adds a degree-aware de-circularized estimator: a fixed-size,
probability-proportional-to-degree (Pareto) null replaces the uniform
permutation, with a bias-corrected enrichment ratio, negative-control empirical
calibration, a conditional portability metric, and a synthetic-graph simulation
that verifies type I error control and power.
"""

__version__ = "2.1.0"

__all__ = ["__version__"]
