%% UC-PBM v6.0

% Developed and written by Firnaaz Ahamed (2020)
% Contact: firnaaz.ahamed@monash.edu; firnaaz.ahamed@gmail.com

%% Simulation Code

close all; clear all; clc;

% Setting input parameters

% Cellulose, cellobiose, biomass loadings

mss = 4.59;                 % Initial concentration of cellulose, g/L
cellobiose = 1e-4;          % Initial concentration of cellobiose, g/L
biomass = 0.04;             % Initial concentration of biomass, g/L

Av_pc = 10;                 % Pre-culture Avicel concentration, g/L
CB_pc = 0;                  % Pre-culture Cellobiose concentration, g/L

% Cellulosome, cellulase, biomass settings (if biomass grown on cellulose)

MW_cellome = 1.73e6;        % Molecular weight of cellulosomes, g/mol
mass_frac_cellulase = 0.50; % Fraction of cellulase in cellulosome
mol_exo = 11.08*0.99;       % Moles of exo-enzymes per mole of cellulosome
mol_endo = 11.08-mol_exo;   % Moles of endo-enzymes per mol of cellulosome
biom_DCW_prtn = 1/1.2;      % Fraction of biomass in DW (N20%)

% Cellulosome, cellulase, biomass settings (if biomass grown on cellobiose)

% MW_cellome = 1.89e6;
% mass_frac_cellulase = 0.39; 
% mol_exo = 9.48*0.99;
% mol_endo = 9.48-mol_exo;
% biom_DCW_prtn = 1/1.02;   % (N2%)

% Mesh setting

p = 20;                    % Number of discrete pivots 

% Time Span 

t_end = 1.62e5;             % Time span of simulation, seconds

% Initial Distribution - base fitted cellulose distribution

load('Raw_Avicel_calibrated_initialDist_Engel_2012_final','MnL','MwL',...
    'MnH','MwH','N','r_massL');

% EFVs from metabolic network

load EFM_NewReduced9aa_CTh_DSM1313_iAT601_Thompson_2016.mat;    % (N20%)
% load EFM_NewReduced10_CTh_DSM1313_iAT601_Thompson_2016.mat;   % (N2%)

% List of Model Paramaters

% ML-PBM parameters

% Exo-enzyme-related parameters

% 1  kp_h_exo - Rate contant of exo-enzyme hydrolysis, 1/s
% 2  kp_f_exo - Rate constant of exo-enzyme complexation, DP.L/mol.s
% 3  kp_e_exo - Rate constant of exo-enzyme decomplexation, 1/s
% 4  tau_exo - Exo-enzyme footprint/area coverage, mol/m2
% 5  k_ads_exo - Rate constant of exo-enzyme adsorption, L/mol.s
% 6  k_des_exo - Rate constant of exo-enzyme desoprtion, 1/s
% 7  k_If_exo(1) - Forward glucose inhibition rate constant (exo-enzyme), L/mol.s
% 8  k_If_exo(2) - Forward cellobiose inhibition rate constant(exo-enzyme), L/mol.s
% 9  k_Ir_exo(1) - Reverse glucose inhibition rate constant(exo-enzyme), 1/s
% 10 k_Ir_exo(2) - Reverse cellobiose inhibition rate constant(exo-enzyme), 1/s

% Substrate-related parameters

% 11 pnt_zone - Number of layers in penetration zone
% 12 r_massL - Mass ratio of cellulose in penetration zone

% Endo-enzyme-related parameters

% 13 kp_h_endo - Rate constant of endo-enzyme hydrolysis (insoluble cellulose), DP/s 
% 14 kp_hs_endo - Rate constant of endo-enzyme hydrolysis (soluble cellulose), L/mol.s
% 15 kp_f_endo - Rate constant of endo-enzyme complexation, DP.L/mol.s
% 16 kp_e_endo - Rate constant of endo-enzyme decomplexation, 1/s
% 17 tau_endo - Endo-enzyme enzyme footprint/area coverage, mol/m2
% 18 k_ads_endo - Rate constant of endo-enzyme adsorption, L/mol.s
% 19 k_des_endo - Rate constant of endo-enzyme desorption, 1/s
% 20 k_If_endo(1) - Forward glucose inhibition rate constant (endo-enzyme), L/mol.s
% 21 k_If_endo(2) - Forward cellobiose inhibition rate constant (endo-enzyme), L/mol.s
% 22 k_Ir_endo(1) - Reverse glucose inhibition rate constant (endo-enzyme), 1/s
% 23 k_Ir_endo(2) - Reverse cellobiose inhibition rate constant (endo-enzyme), 1/s

% L-HCM parameters

% 24 aF(1) - Constitutive intracellular enzyme synthesis rate (F1), 1/s
% 25 aF(2) - Constitutive intracellular enzyme synthesis rate (F2), 1/s
% 26 bF(1) - Intracellular enzyme decay rate (F1),1/s
% 27 bF(2) - Intracellular enzyme decay rate (F2),1/s
% 28 kE(1) - Inducive intracellular enzyme synthesis rate (F1), 1/s
% 29 kE(2) - Inducive intracellular enzyme synthesis rate (F2), 1/s
% 30 KE(1) - MM constant for intracellular enzyme synthesis (F1), mol/L
% 31 KE(2) - MM constant for intracellular enzyme synthesis (F2), mol/L
% 32 k_inh_G - Glucose inhibition constant, mol/L
% 33 k_inh_etoh - Ethanol inhibition constant, mol/L
% 34 kmax(1) - Maximum substrate uptake rate (F1), mol/g-biom.s
% 35 kmax(2) - Maximum substrate uptake rate (F2), mol/g-biom.s
% 36 K(1) - MM constant for substrate uptake (F1), mol/L
% 37 K(2) - MM constant for substrate uptake (F2), mol/L
% 38 aE - Constitutive cellulase synthesis rate, g/g-biom.s

% Input the values for model parameters below:

prm(1) = 2.15;
prm(2) = 3e7;
prm(3) = 10;
prm(4) = 2.22e-8;
prm(5) = 7.3722e5;
prm(6) = 0.01;
prm(7) = 0;
prm(8) = 0;
prm(9) = 0;
prm(10) = 0;

prm(11) = 1;
prm(12) = 0.028;

prm(13) = 4;
prm(14) = 3;
prm(15) = 1e9;
prm(16) = 10;
prm(17) = 8.04e-9;
prm(18) = 7.11e5;
prm(19) = 0.01;
prm(20) = 0;
prm(21) = 0;
prm(22) = 0;
prm(23) = 0;

prm(24) = 0.1/3600;
prm(25) = 0.1/3600;
prm(26) = 0.2/3600;
prm(27) = 0.2/3600;
prm(28) = 1/3600;	
prm(29) = 1/3600;	
prm(32) = 1.22e-3;	
prm(33) = 1.74;     
prm(34) = 1.3e-6;	
prm(35) = 1.9e-6;	
prm(36) = 1e-6;	
prm(37) = 9e-6;	
prm(30) = prm(36);	
prm(31) = prm(37);	
prm(38) = 3e-7;

% Run simulation

[Z,idx_cb,idx_atp,idx_csm,idx_biom,idx_etoh,idx_lac,idx_form,idx_ace,...
    F1,F1A,F1C,F2,F2A,F2B,idx_met,ZF1,ZF2,Y1_met_lump,Y2_met_lump,Y1_E,...
    Y2_E,Y2_B,x_piv,R_in,layers,C_T,C_S,C_INT,MnT,MwT,pdiT,t,C,CXS,CNS,...
    CNXS,CFS,CS,CI,CT,CC,E_T_exo,E_T_endo,CB,etoh,lac,form,ace,cellome,...
    biom,cellulase_tot,cellu,cell_gluc_eq,eF_rel,R,conv_rem,conv_cb,uF,vF,...
    r_up,mu] = UCPBM(p,N,prm,mss,MnL,MwL,MnH,MwH,t_end,CB_pc,cellobiose,...
    biomass,MW_cellome,mol_exo,mol_endo,biom_DCW_prtn,ems,ems_idx,reac,...
    mass_frac_cellulase);

% List of outputs

% x_piv     % DP mesh
% R_in      % Initial particle radius
% layers    % Total no. of layers in cellulose particles
% C_T       % Total initial cellulose distribution
% C_S       % Surface initial cellulose distribution
% C_INT     % Internal initial cellulose distribution
% MnT       % Overall initial number-average DP
% MwT       % Overall initial weight-average DP
% pdiT      % Overall initial polydispersity index

% Indexes of key metabolic rxns

% idx_cb,idx_atp,idx_csm,idx_biom,idx_etoh,idx_lac,idx_form,idx_ace

% Metabolic information

% Z                                         % EFV matrix
% F1,F1A,F1C,F2,F2A,F2B                     % Number of EFVs in EFV families
% ZF1,ZF2                                   % Lumped EFV families
% Y1_met_lump,Y2_met_lump,Y1_E,Y2_E,Y2_B    % Lumped yields

% t         % Time, s

% Time,DP-dependent variables

% C         % Soluble products, mol/L
% CXS       % CBH-bound complex, mol/L
% CNS       % EG-bound complex, mol/L
% CNXS      % CBH-EG-bound complex, mol/L
% CFS       % Free un-bound surface polymers, mol/L
% CS        % Total surface polymers, mol/L
% CI        % Internal polymers, mol/L 
% CT        % Total polymers, mol/L
% CC        % Insoluble polymers, mol/L

% Time-dependent variables

% E_T_endo      	% Total endo-enzymes, mol/L
% E_T_exo           % Total exo-enzymes, mol/L
% CB                % Cellobiose, g/L
% etoh              % Ethanol, g/L
% lac               % Lactate, g/L
% form              % Formate, g/L
% ace               % Acetate, g/L
% cellome           % Cellulosome, g/L
% biom              % Biomass, g/L
% cellulase_tot     % Total cellulase, g/L 
% cellu             % Remaining insoluble cellulose, g/L
% cell_gluc_eq      % Remaining insoluble cellulose, g glu/L
% eF_rel            % Relative intracellular enzyme levels
% R                 % Transient of particle radius, m
% conv_rem          % Overall cellulose conversion
% conv_cb           % Overall cellobiose conversion
% uF,vF             % Cybernetic variables
% r_up              % Uptake fluxes through EFV families, mol/g-biom.s
% mu                % Specific biomass growth rate, 1/s