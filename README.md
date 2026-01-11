# Project Overview

This repository contains the full solution for the MRI Scheduling Case Study of the course Computational Research Skills. The project consists of two main components: 

* Statistical analysis of historical MRI data (Part 1).

* Discrete-event simulation and policy evaluation (Part 2).

The objective is to translate empirical data into reliable inputs and assess operational performance under different scheduling strategies. 

# Case Description 

* Type 1 patients: The number of arrivals (calls to schedule a MRI appointment) follows a Poisson distribution with unknown mean.
* Type 2 patients: No parametric distribution is assumed for arrivals or scan durations.
* The hospital is primarily interested in interpretable operational quantities, rather than model parameters alone.
* Uncertainty in estimates is quantified using bootstrap methods, and robustness is assessed through simulation.
Only working hours (08:00-17:00) are considered. Time outside working horus is treated as not passing.

# Data
Input file: `ScanRecords.csv`. 
The file contains historical MRI scan request data with the following variables: 
* `Date`: Date of the scan request.
* `Time`: Time of the request (in hours).
* `Duration`: Scan duration (in hours).
* `PatientType`: Patient category (`Type 1` or `Type 2`).

# Part 1: Statistical Analysis
The R script(s) performs: 
1. Data processing
* Cleaning and validation of observations.
* Construction of timestamps and working-time variables.
* Seperation by patient type. 
2. Exploratory analysis
* Daily arrival patterns.
* Interarrival times during working hours.
* Scan duration distributions.
3. Estimation of operational quantities
* Mean arrivals per day.
* Quantiles of arrivals and durations.
* Probabilities of long scan durations.
4. Uncertainty quantification
* Parametric boostrap for Type 1 arrivals and durations.
* Nonparametric boostrap for Type 2 arrivals and durations.
5. Robustness analysis (Type 2)
* Monte Carlo simulation under alternative data-generating models.
* Evaluation of coverage, bias, and root mean squared error.
6. Export of simulation input
* Estimated paramters and empirical distributions are exported as CSV files for use in the simulation model.

All statistical results and interpretations are reported in the report handed in via Canvas. 

# Part 2: Simulation and Policy Evaluation
The simulation component implements a discrete-event simulation model of the MRI scheduling system. It uses the CSV files generated in Part 1 as inputs for:
* Arrival processes of both patient types.
* Scan duration distributions.
* Working-hour constraints and slot lengths.
Multiple scheduling policies are evaluated and compared based on operational performance metrics.

# Authors
* Alp Erçikan
* Qiyuan Li
* Laura Maas
* Mart van der Vleuten
