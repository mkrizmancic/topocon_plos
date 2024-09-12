# Distributed topology control of multi-agent systems

This repository contains the code for the paper "_Adaptive connectivity control in networked multi-agent systems: A distributed approach_" by M. Krizmancic and S. Bogdan submitted to PLOS ONE journal on April 30, 2024.

Also check out the final dataset and figures generated with this code: https://zenodo.org/doi/10.5281/zenodo.11082641

## Abstract
Effective communication is crucial for the performance and collaboration within cooperative networked multi-agent systems. However, existing literature lacks comprehensive solutions for dynamically monitoring and adjusting communication topologies to balance connectivity and energy efficiency. This study addresses this gap by proposing a distributed approach for estimating and controlling system connectivity over time.
We introduce a modified consensus protocol where agents exchange local assessments of communication link quality, enabling the estimation of a global weighted adjacency matrix without requiring centralized information. The system's connectivity is measured using the second smallest eigenvalue of the communication graph Laplacian, commonly referred to as algebraic connectivity. Additionally, we enhance the consensus protocol with an adaptive mechanism to expedite convergence, irrespective of system size or structure.
Furthermore, we present an analytical method for connectivity control based on the Fiedler vector approximation, facilitating the addition or removal of communication links. This method adjusts control parameters to accommodate minor variations in link quality while reconfiguring the network in response to significant changes. Notably, it identifies and eliminates energy-consuming yet non-contributory links, improving long-term connectivity efficiency.
Simulation experiments across diverse scenarios and the number of agents validate the efficacy of our proposed algebraic connectivity estimation and tracking strategy. Results demonstrate robust connectivity maintenance against external disturbances and agent failures, underscoring the practical utility of our approach for real-world multi-agent systems.

## Code organization
The code was written in MATLAB R2023b.

| **Directory** | **Description**                                                     |
|---------------|---------------------------------------------------------------------|
| helper        | General purpose helper functions                                    |
| plotting      | Functions used for plotting                                         |
| procedures    | Functions for essential steps of the algorithm                      |
| results       | Directory for storing figures and mat files                         |
| runtime       | Functions for simulating external effects during runtime            |
| setup         | Loading and generating initial topologies and simulation parameters |
| test          | Functions for testing various approaches and methods                |


## Usage

To replicate the figures in the paper:
1. Set the desired figure number in line 13 of the `main.m` script.
2. Run the `main.m` script.