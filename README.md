stochastic
==========

This is my personal PhD research code.

It performs 1D neutron radiation transport calculations.

It can perform classical Monte Carlo calculations in discrete material chunks based upon random chord sampling, Levermore-Pomraning methods, or an atomic mix of input.

It can also compute transport using Woodcock Monte Carlo on continuous material cross sections generated from various approaches to the Karhunen-Loeve expansion.

Included are KL tools, stochastic collocation functionality, and soon to be polynomial chaos functionality.

At various times chronicled in commits, it also incorporated a basic diffusion and Sn deterministic code, performed Multi-level Monte Carlo, and performed what I have called "Weight-adjusted Monte Carlo", for which distance to collision is sampled using a fictitious total cross section (enabling even negative cross sections) and particle weights are adjusted to maintain statistics.
