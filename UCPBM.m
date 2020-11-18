function [Z,idx_cb,idx_atp,idx_csm,idx_biom,idx_etoh,idx_lac,idx_form,...
    idx_ace,F1,F1A,F1C,F2,F2A,F2B,idx_met,ZF1,ZF2,Y1_met_lump,...
    Y2_met_lump,Y1_E,Y2_E,Y2_B,x_piv,R_in,layers,C_T,C_S,C_INT,MnT,MwT,...
    pdiT,t,C,CXS,CNS,CNXS,CFS,CS,CI,CT,CC,E_T_exo,E_T_endo,CB,etoh,lac,...
    form,ace,cellome,biom,cellulase_tot,cellu,cell_gluc_eq,eF_rel,R,...
    conv_rem,conv_cb,uF,vF,r_up,mu] = UCPBM(p,N,prm,mss,MnL,MwL,MnH,MwH,...
    t_end,CB_pc,cellobiose,biomass,MW_cellome,mol_exo,mol_endo,...
    biom_DCW_prtn,ems,ems_idx,reac,mass_frac_cellulase)

[Z,idx_cb,idx_atp,idx_csm,idx_biom,idx_etoh,idx_lac,idx_form,idx_ace,F1,...
    F1A,F1C,F2,F2A,F2B,idx_met,ZF1,ZF2,Y1_met_lump,Y2_met_lump,Y1_E,...
    Y2_E,Y2_B] = EFV_setting(ems,ems_idx,biom_DCW_prtn,reac);

[x_piv,R_in,layers,C_T,C_S,C_INT,MnT,MwT,pdiT,q,vm,R0,ro,L,n,exp_ratio]...
    = MLPBM_setting(p,N,prm,mss,MnL,MwL,MnH,MwH);

[t,y,eF_max,kmax,K,K_inh_G,K_inh_etoh] = ODE_setting(p,q,C_S,C_INT,prm,...
    x_piv,vm,t_end,Y1_met_lump,Y2_met_lump,Y1_E,Y2_E,Y2_B,CB_pc,...
    cellobiose,biomass,MW_cellome,mol_exo,mol_endo,R0,ro,L,n,exp_ratio);

[C,CXS,CNS,CNXS,CFS,CS,CI,CT,CC,E_T_exo,E_T_endo,CB,etoh,lac,form,ace,...
    cellome,biom,cellulase_tot,cellu,cell_gluc_eq,eF_rel,R,conv_rem,...
    conv_cb,uF,vF,r_up,mu] = post_processing(t,y,p,q,MW_cellome,mol_endo,...
    mol_exo,mass_frac_cellulase,x_piv,eF_max,ro,L,n,mss,cellobiose,kmax,...
    K,K_inh_G,K_inh_etoh,Y2_B);

end

%%-------------------------------------------------------------------------

function [Z,idx_cb,idx_atp,idx_csm,idx_biom,idx_etoh,idx_lac,idx_form,...
    idx_ace,F1,F1A,F1C,F2,F2A,F2B,idx_met,ZF1,ZF2,Y1_met_lump,...
    Y2_met_lump,Y1_E,Y2_E,Y2_B] = EFV_setting(ems,ems_idx,biom_DCW_prtn,...
    reac)
    
% EFV setting

n_eta = 1;                      % Sensitivity setting to tuning parameters

% Tuning parameters for N20%

a_tune1 = [-1.544065229170075 1.569902586838220...
    4.730079923633626 -8.058648461794135];
a_tune2 = [-0.002019303084171 -0.007859218375312...
    -0.009091575956450 0.005826206931968];

% Tuning parameters for N2%

% a_tune1 = [-4.467191370597475 4.577593783957661...
%     13.788127542564130 -23.705397747646195];        
% a_tune2 = [3.644003789716016 -2.686629088044898...
%     -9.470963287592038 13.097851990509195];

Z = ems';                                       % EFV matrix, Z                   
rxn = reac(ems_idx);                            % List of metabolic rxns

% Indexes for key metabolic rxns

idx_cb = find(~cellfun(@isempty,strfind(rxn,'T_e_to_c_C00185_c')));
idx_atp = find(~cellfun(@isempty,strfind(rxn,'R_R_MAINT')));
idx_biom = find(~cellfun(@isempty,strfind(rxn,'T_c_to_e_m85')));
idx_csm = find(~cellfun(@isempty,strfind(rxn,'R_EXC_OUT_m90')));
idx_etoh = find(~cellfun(@isempty,strfind(rxn,'T_c_to_e_C00469_c')));
idx_lac = find(~cellfun(@isempty,strfind(rxn,'T_c_to_e_C00186_c')));
idx_form = find(~cellfun(@isempty,strfind(rxn,'T_c_to_e_C00058_c')));
idx_ace = find(~cellfun(@isempty,strfind(rxn,'T_c_to_e_C00033_c')));

% Classification of EFV families

F = find(Z(idx_cb,:)~=0)';            
F1 = intersect(F,setdiff(find(Z(idx_csm,:)~=0),find(Z(idx_biom,:)~=0)));
F2 = setdiff(F,F1);

e1C = intersect(F1,find(Z(idx_csm,:)~=0));
e1A = intersect(F1,find(Z(idx_atp,:)~=0));  

e2B = intersect(F2,find(Z(idx_biom,:)~=0)); 
e2A = intersect(F2,find(Z(idx_atp,:)~=0));  

F1C = e1C;                                  
F1A = setdiff(e1A,e1C);                     
F1 = unique([F1C;F1A]);                     

F2B = e2B;                                  
F2A = setdiff(e2A,e2B);
F2 = unique([F2B;F2A]);

z1C = Z(:,F1C);                 % z-matrix of F1
z2B = Z(:,F2B);                 % z-matrix of F2,biom    
z2A = Z(:,F2A);                 % z-matrix of F2,atp

% Yields of metabolites in individual EFVs before lumping

Y1C_csm = (z1C(idx_csm,:)./z1C(idx_cb,:));
Y1C_etoh = z1C(idx_etoh,:)./z1C(idx_cb,:);
Y1C_lac = z1C(idx_lac,:)./z1C(idx_cb,:);
Y1C_form = z1C(idx_form,:)./z1C(idx_cb,:);
Y1C_ace = z1C(idx_ace,:)./z1C(idx_cb,:);

Y2B_biom = (z2B(idx_biom,:)./z2B(idx_cb,:));
Y2B_etoh = z2B(idx_etoh,:)./z2B(idx_cb,:);
Y2B_lac = z2B(idx_lac,:)./z2B(idx_cb,:);
Y2B_form = z2B(idx_form,:)./z2B(idx_cb,:);
Y2B_ace = z2B(idx_ace,:)./z2B(idx_cb,:);

Y2A_atp = (z2A(idx_atp,:)./z2A(idx_cb,:));
Y2A_etoh = z2A(idx_etoh,:)./z2A(idx_cb,:);
Y2A_lac = z2A(idx_lac,:)./z2A(idx_cb,:);
Y2A_form = z2A(idx_form,:)./z2A(idx_cb,:);
Y2A_ace = z2A(idx_ace,:)./z2A(idx_cb,:);

% Metabolite order: eth,lac,form,ace

idx_met = [idx_etoh;idx_lac;idx_form;idx_ace];  % Index of metabolite secreting rxns

Y1C_met = [Y1C_etoh;Y1C_lac;Y1C_form;Y1C_ace];
Y2B_met = [Y2B_etoh;Y2B_lac;Y2B_form;Y2B_ace];
Y2A_met = [Y2A_etoh;Y2A_lac;Y2A_form;Y2A_ace];

% EFV lumping scheme

wt = 1;     % Arbitrary weightage

eta1C = zeros(length(F1C),1);

for i = 1:length(F1C)
    eta1C(i) = wt*(Y1C_csm(i)+sum(a_tune1'.*Y1C_met(:,i)));
end

eta2B = zeros(length(F2B),1);

for i = 1:length(F2B)
    eta2B(i) = wt*(Y2B_biom(i)+sum(a_tune2'.*Y2B_met(:,i)));
end

eta2A = zeros(length(F2A),1);

for i = 1:length(F2A)
    eta2A(i) = wt*(Y2A_atp(i)+sum(a_tune2'.*Y2A_met(:,i)));
end

% Lumped EFV matrices

ZF1C = zeros(length(rxn),1);
ZF2B = zeros(length(rxn),1); 
ZF2A = zeros(length(rxn),1); 

for i = 1:length(rxn)
    ZF1C(i) = sum(z1C(i,:)'.*(eta1C.^(3*n_eta)))/sum(eta1C.^(3*n_eta));
    ZF2B(i) = sum(z2B(i,:)'.*(eta2B.^(3*n_eta)))/sum(eta2B.^(3*n_eta));
    ZF2A(i) = sum(z2A(i,:)'.*(eta2A.^(3*n_eta)))/sum(eta2A.^(3*n_eta));
end

ZF1 = ZF1C;
ZF2 = ZF2B+ZF2A;

% Final yields of metabolites after EFV lumping

Y1_met_lump = zeros(length(idx_met),1);                 % mol/mol
Y2_met_lump = zeros(length(idx_met),1);

for i = 1:length(idx_met)
    Y1_met_lump(i) = (ZF1(idx_met(i))/ZF1(idx_cb));
    Y2_met_lump(i) = (ZF2(idx_met(i))/ZF2(idx_cb));
end

Y1_E = (ZF1(idx_csm)/ZF1(idx_cb))*1000;                 % g/mol
Y2_E = (ZF2(idx_csm)/ZF2(idx_cb))*1000;

Y2_B = biom_DCW_prtn*(ZF2(idx_biom)/ZF2(idx_cb))*1000;  % g/mol

end

%%-------------------------------------------------------------------------

function [x_piv,R_in,layers,C_T,C_S,C_INT,MnT,MwT,pdiT,q,vm,R0,ro,L,n,...
    exp_ratio] = MLPBM_setting(p,N,prm,mss,MnL,MwL,MnH,MwH)

% General settings

vm = 1;   % Size of monomer

qmax = 1+(log(N/(p+1))/log(1+(vm/(p+1))));  % Maximum pivots in cont region
q = floor(qmax);                            % No. of pivots in cont region
ratio = (N/(p+1))^(1/(q-1));                % Ratio of geometric progr

% Pivots for discrete-continuous mesh

x_piv = zeros(p+q,1);
x_piv(1:p) = 1:p;

for i = p+1:p+q
    x_piv(i) = (p+1)*ratio^(i-(p+1));       % Geometric mesh
end

% Boundary points for discrete-continuous mesh

x_bound = zeros(p+q+1,1);
x_bound(1) = 0.5;

for i=2:p+q
    x_bound(i) = (x_piv(i)+x_piv(i-1))/2;
end

x_bound(p+q+1) = x_piv(end);

% Initial distribution

c_in = @(x,alpha,beta,pin) pin*gampdf(x,alpha,beta); % Gamma distribution

% Physical properties

layers = round(2/(((mss*prm(12))/prm(11))/mss));     % Total no. of layers
ro = 1500;                   % Cellulose density, g/L
L = (N/2)*10.38e-10;         % Length of microfibril, m
R0 = 1e-9;                   % Width of single layer, m
R_in = layers*R0;            % Total radius of particle, m
n = mss/(ro*pi*(R_in^2)*L);  % No. of particles per unit vol, 1/m3
p_zone = prm(11);            % No. of layers in penetration zone
i_zone = (R_in/R0)-p_zone;   % No. of layers in internal zone

% Radiuses of discrete layers

R_set = sort(0:R0:R_in,'descend')';

% Mass of polymers in each layer

mass_l = zeros(length(R_set)-1,1); 

for i = 2:length(R_set)
    mass_l(i-1) = (n*ro*pi*(R_set(i-1)^2)*L)-(n*ro*pi*(R_set(i)^2)*L);
end

% Initial cellulose distribution in each layer

pin_l = zeros(length(mass_l),1);
c_in_tl = zeros(length(mass_l),p+q); % Molar concentration density
C_in_tl = zeros(length(mass_l),p+q); % Molar concentration

optionsfmincon = optimoptions('fsolve','display','none');

for i = 1:p_zone
    pin_l(i) = fsolve(@(x) sum(c_in(x_piv,MnL/(MwL-MnL),MwL-MnL,x)...
        .*(x_bound(2:end)-x_bound(1:p+q)).*(162*x_piv+18))-mass_l(i),0,...
        optionsfmincon);
    c_in_tl(i,:) = c_in(x_piv,MnL/(MwL-MnL),MwL-MnL,pin_l(i));
end

for i = p_zone+1:p_zone+i_zone
    pin_l(i) = fsolve(@(x) sum(c_in(x_piv,MnH/(MwH-MnH),MwH-MnH,x)...
        .*(x_bound(2:end)-x_bound(1:p+q)).*(162*x_piv+18))-mass_l(i),0,...
        optionsfmincon);
    c_in_tl(i,:) = c_in(x_piv,MnH/(MwH-MnH),MwH-MnH,pin_l(i));
end

for i = 1:p+q
    C_in_tl(:,i) = c_in_tl(:,i).*(x_bound(i+1)-x_bound(i));
end

% Overall initial cellulose distribution

C_T = zeros(p+q,1);

for i = 1:p+q
    C_T(i) = sum(C_in_tl(:,i));
end

% Surface initial cellulose distribution

C_S = zeros(p+q,1);

for i = 1:p+q
    C_S(i) = sum(C_in_tl(1,i));
end

% Internal initial cellulose distribution

C_INT = zeros(p+q,1);

for i = 1:p+q
    C_INT(i) = sum(C_in_tl(2:end,i));
end

% Ratio of molar concentration to total mass in each layer 

C_ratio = zeros(length(R_set)-1,p+q);

for i = 1:length(R_set)-1
    for j = 1:p+q
        C_ratio(i,j) = C_in_tl(i,j)/sum(C_in_tl(i,:)'.*(162*x_piv+18));
    end
end

% Piece-wise function for ratio of molar concentration to total mass as a
% function of particle radius

exp_ratio = cell(p+q,1);

for i = 1:p+q
    exp_ratio{i} = mkpp(flipud(R_set),flipud(C_ratio(:,i)));
end

% General properties of initial cellulose distribution

MnT = sum(C_T.*x_piv)/sum(C_T);
MwT = sum(C_T.*(x_piv.^2))/sum(C_T.*x_piv);
pdiT = MwT/MnT;

end

%%-------------------------------------------------------------------------

function [t,y,eF_max,kmax,K,K_inh_G,K_inh_etoh] = ODE_setting(p,q,C_S,...
    C_INT,prm,x_piv,vm,t_end,Y1_met_lump,Y2_met_lump,Y1_E,Y2_E,Y2_B,...
    CB_pc,cellobiose,biomass,MW_cellome,mol_exo,mol_endo,R0,ro,L,n,...
    exp_ratio)

% Parameter setting

% ML-PBM parameters

% Exo-enzyme-related parameters

kp_h_exo = prm(1);                      % Rate constant of hydrolysis, 1/s     
m_h_exo = 0;
k_h_exo = kp_h_exo*(x_piv.^m_h_exo);    % Rate coefficients of hydrolysis, 1/s

kp_f_exo = prm(2);                      % Rate constant of complexation, DP.L/mol.s
m_f_exo = -1;
k_f_exo = kp_f_exo*(x_piv.^m_f_exo);    % Rate coefficients of complexation, L/mol.s

kp_e_exo = prm(3);                      % Rate constant of decomplexation, 1/s
m_e_exo = 0;
k_e_exo = kp_e_exo*(x_piv.^m_e_exo);    % Rate coefficients of decomplexation, 1/s

tau_exo = prm(4);              % Enzyme footprint, mol/m2  
k_ads_exo = prm(5);            % Adsorption constant, L/mol.s
k_des_exo = prm(6);            % Desorption constant, 1/s

k_If_exo = [prm(7) prm(8)];    % Forward glucose/cellobiose inhibition, L/mol.s 
k_Ir_exo = [prm(9) prm(10)];   % Reverse glucose/cellobiose inhibition, 1/s 

% Endo-enzyme-related parameters

kp_h_endo = prm(13);                    % Rate constant of hydrolysis (insoluble cellulose), 1/s
m_h_endo = 0;                           % Endo-enzymes in C.th unaffected by crystallinity
k_h_endo = kp_h_endo*(x_piv.^m_h_endo); % Rate coefficients of hydrolysis(insoluble cellulose), 1/s

kp_hs_endo = prm(14);           % Rate constant of hydrolysis (soluble cellulose), 1/s         
m_hs_endo = 0;
k_hs_endo = kp_hs_endo*(x_piv.^m_hs_endo);  
k_hs_endo(2) = 0;               % Rate coefficients of hydrolysis (soluble cellulose), L/mol.s
                                % Endo-enzymes in C.th do not hydrolyze cellobiose

kp_f_endo = prm(15);                     % Rate constant of complexation, L/mol.s   
m_f_endo = 0;                            % Endo-enzymes in C.th unaffected by crystallinity
k_f_endo = kp_f_endo*(x_piv.^m_f_endo);  % Rate coefficients of complexation, L/mol.s

kp_e_endo = prm(16);                     % Rate constant of decomplexation, 1/s
m_e_endo = 0;
k_e_endo = kp_e_endo*(x_piv.^m_e_endo);  % Rate coefficients of decomplexation, 1/s

tau_endo = prm(17);                      % Enzyme footprint, mol/m2 
k_ads_endo = prm(18);                    % Adsorption constant, L/mol.s  
k_des_endo = prm(19);                    % Desorption constant, 1/s

k_If_endo = [prm(20) prm(21)];           % Forward glucose/cellobiose inhibition, L/mol.s 
k_Ir_endo = [prm(22) prm(23)];           % Reverse glucose/cellobiose inhibition, 1/s 

% L-HCM parameters                

aF(1) = prm(24);        % Constitutive intracellular enzyme synthesis rate (F1), 1/s
aF(2) = prm(25);        % Constitutive intracellular enzyme synthesis rate (F2), 1/s
bF(1) = prm(26);        % Intracellular enzyme decay rate (F1),1/s
bF(2) = prm(27);        % Intracellular enzyme decay rate (F2),1/s
kE(1) = prm(28);        % Inducive intracellular enzyme synthesis rate (F1), 1/s
kE(2) = prm(29);        % Inducive intracellular enzyme synthesis rate (F2), 1/s
KE(1) = prm(30);        % MM constant for intracellular enzyme synthesis (F1), mol/L
KE(2) = prm(31);        % MM constant for intracellular enzyme synthesis (F2), mol/L
K_inh_G = prm(32);      % Glucose inhibition constant, mol/L
K_inh_etoh = prm(33);   % Ethanol inhibition constant, mol/L
kmax(1) = prm(34);      % Maximum substrate uptake rate (F1), mol/g-biom.s
kmax(2) = prm(35);      % Maximum substrate uptake rate (F2), mol/g-biom.s
K(1) = prm(36);         % MM constant for substrate uptake (F1), mol/L
K(2) = prm(37);         % MM constant for substrate uptake (F2), mol/L
aE = prm(38);           % Constitutive cellulase synthesis rate, g/g-biom.s

% Maximum intracellular enzyme levels

eF_max(1) = (aF(1)+kE(1))/bF(1);
eF_max(2) = (aF(2)+kE(2))/(bF(2)+(Y2_B*kmax(2)));

% Initial intracellular enzyme levels

[init_e] = int_enz(aF,kE,bF,Y2_B,kmax,CB_pc,K,KE);

% Fixed pivot technique

[n_exo,n_endo] = fp(p,q,vm,x_piv);      % Particle allocation functions

% Initial conditions

init = zeros(6*p+6*q+17,1);
init(1:p+q) = 0;                        % 1 to p+q --> free liberated polymers
init(2) = cellobiose/342;               % 2 --> cellobiose
init(p+q+1:2*p+2*q) = 0;                % p+q+1 to 2p+2q --> CBH-bound polymers
init(2*p+2*q+1:3*p+3*q) = 0;            % 2p+2q+1 to 3p+3q --> EG-bound polymers
init(3*p+3*q+1:4*p+4*q) = 0;            % 3p+3q+1 to 4p+4q --> EG-CBH-bound polymers
init(4*p+4*q+1:5*p+5*q) = C_S;          % 4p+4q+1 to 5p+5q --> surface-accessible polymer
init(5*p+5*q+1:6*p+6*q) = C_INT;        % 5p+5q+1 to 6p+6q --> internal inaccessible polymer
init(6*p+6*q+1) = 0;                    % 6p+6q+1 --> surface adsorbed EG
init(6*p+6*q+2) = 0;                    % 6p+6q+2 --> surface adsorbed CBH
init(6*p+6*q+3:6*p+6*q+4) = 0;          % 6p+6q+3 to 6p+6q+4 --> inhibited EG
init(6*p+6*q+5:6*p+6*q+6) = 0;          % 6p+6q+5 to 6p+6q+6 --> inhibited free CBH
init(6*p+6*q+7:6*p+6*q+8) = 0;          % 6p+6q+7 to 6p+6q+8 --> inhibited adsorbed CBH    
init(6*p+6*q+9) = 0;                    % 6p+6q+9 --> inhibited BG (placeholder for future extension)
init(6*p+6*q+10) = init_e(1)*eF_max(1); % 6p+6q+10 --> intracellular enzyme (F1)
init(6*p+6*q+11) = init_e(2)*eF_max(2); % 6p+6q+11 --> intracellular enzyme (F2)
init(6*p+6*q+12) = 0;                   % 6p+6q+12 --> etoh
init(6*p+6*q+13) = 0;                   % 6p+6q+13 --> lac
init(6*p+6*q+14) = 0;                   % 6p+6q+14 --> form
init(6*p+6*q+15) = 0;                   % 6p+6q+15 --> ace
init(6*p+6*q+16) = 0;                   % 6p+6q+16 --> cellome
init(6*p+6*q+17) = biomass;             % 6p+6q+17 --> biomass

% ODE solver

tstart = 0;                                     
tfinal = t_end;                   

tic;

options = odeset('AbsTol',1e-10,'RelTol',1e-10,'NonNegative',1:(6*p+6*q+17),...
    'OutputFcn',@odeprog,'Events',@odeabort);

[t,y] = ode15s(@(t,y) UCPBM_ode(p,q,n_exo,n_endo,k_h_exo,k_h_endo,k_hs_endo,...
    k_f_exo,k_f_endo,k_e_exo,k_e_endo,k_ads_endo,k_des_endo,k_ads_exo,...
    k_des_exo,tau_endo,tau_exo,k_If_exo,k_Ir_exo,k_If_endo,k_Ir_endo,...
    x_piv,R0,ro,L,n,exp_ratio,kmax,K,kE,K_inh_G,K_inh_etoh,eF_max,...
    Y2_B,aF,bF,Y1_met_lump,Y2_met_lump,aE,Y1_E,...
    Y2_E,MW_cellome,mol_exo,mol_endo,y,t),[tstart tfinal],init,options);

time_fp_sim = '\nTotal time elapsed for UC-PBM simulation is %.1f seconds.\n';

fprintf(time_fp_sim,toc)

end

%%-------------------------------------------------------------------------

function [init_e] = int_enz(aF,kE,bF,Y2_B,kmax,CB_pc,K,KE)

e_rel_lb(1) = aF(1)/(aF(1)+kE(1));
e_rel_lb(2) = (aF(2)/(bF(2)+Y2_B*kmax(2)*((CB_pc/342)/(K(2)+(CB_pc/342)))))*...
    ((bF(2)+Y2_B*kmax(2))/(aF(2)+kE(2)));
e_rel_ub(1) = (aF(1)+kE(1)*((CB_pc/342)/(KE(1)+(CB_pc/342))))/(aF(1)+kE(1));
e_rel_ub(2) = ((aF(2)+kE(2)*((CB_pc/342)/(KE(2)+(CB_pc/342))))/...
    (bF(2)+Y2_B*kmax(2)*((CB_pc/342)/(K(2)+(CB_pc/342)))))*...
    ((bF(2)+Y2_B*kmax(2))/(aF(2)+kE(2)));

if CB_pc == 0
    init_e(1) = e_rel_ub(1);
    init_e(2) = e_rel_lb(2);
else
    init_e(1) = e_rel_lb(1);
    init_e(2) = e_rel_ub(2);
end

end

%%-------------------------------------------------------------------------

function [n_exo,n_endo] = fp(p,q,vm,x_piv)

% Particle allocation function - fixed pivot technique

% Chain-end dimer scission

n_exo = zeros(p+q,p+q);
n_exo(1,3) = 1;
n_exo(2,3) = 1; n_exo(2,4) = 2; n_exo(2,5:end) = 1;

for i = 3:p+q
    for j = i:p+q
        if j==i && x_piv(j)-2*vm>x_piv(i-1)
            n_exo(i,j) = ((x_piv(j)-2*vm)-x_piv(i-1))/(x_piv(i)-x_piv(i-1));
        elseif j==i && x_piv(j)-2*vm==x_piv(i-1)
            n_exo(i,j) = 0;
        elseif j~=i && x_piv(j)-2*vm>x_piv(i-1) && x_piv(j)-2*vm<x_piv(i)
            n_exo(i,j) = ((x_piv(j)-2*vm)-x_piv(i-1))/(x_piv(i)-x_piv(i-1));
        elseif j~=i && x_piv(j)-2*vm>x_piv(i) && x_piv(j)-2*vm<x_piv(i+1)
            n_exo(i,j) = (x_piv(i+1)-(x_piv(j)-2*vm))/(x_piv(i+1)-x_piv(i));
        elseif j~=i && x_piv(j)-2*vm==x_piv(i)
            n_exo(i,j) = 1;
        elseif j~=i && x_piv(j)-2*vm==x_piv(i-1)
            n_exo(i,j) = 0;
        else
            n_exo(i,j) = 0;
        end
    end
    n_exo(i,i) = n_exo(i,i)-1;
end

% Random scission

n_endo = zeros(p+q,p+q);

for j = 2:p+q
    n_endo(1,j) = 2/(x_piv(j-1));
end

for i = 2:p
    for j = i+1:p+q
        n_endo(i,i) = -1;
        n_endo(i,j) = 2/(x_piv(j-1));
    end
end

for i = p+1:p+q
    for j = i+1:p+q
        n_endo(i,j) = (x_piv(i+1)-x_piv(i-1))/(x_piv(j-1));
    end
    n_endo(i,i) = -1;
end

end

%%-------------------------------------------------------------------------

function ode = UCPBM_ode(p,q,n_exo,n_endo,k_h_exo,k_h_endo,k_hs_endo,...
    k_f_exo,k_f_endo,k_e_exo,k_e_endo,k_ads_endo,k_des_endo,k_ads_exo,...
    k_des_exo,tau_endo,tau_exo,k_If_exo,k_Ir_exo,k_If_endo,k_Ir_endo,...
    x_piv,R0,ro,L,n,exp_ratio,kmax,K,kE,K_inh_G,K_inh_etoh,eF_max,...
    Y2_B,aF,bF,Y1_met_lump,Y2_met_lump,aE,Y1_E,...
    Y2_E,MW_cellome,mol_exo,mol_endo,y,t)

dydt = zeros(6*p+6*q+17,1);

% Enzymes

E_T_endo = mol_endo*(y(6*p+6*q+16)/MW_cellome);
E_T_exo = mol_exo*(y(6*p+6*q+16)/MW_cellome);

E_F_endo = E_T_endo-y(6*p+6*q+1)-sum(y(2*p+2*q+1:4*p+4*q))...
    -sum(y(6*p+6*q+3:6*p+6*q+4));
E_F_exo = E_T_exo-y(6*p+6*q+2)-sum(y(p+q+1:2*p+2*q))...
    -sum(y(3*p+3*q+1:4*p+4*q))-sum(y(6*p+6*q+5:6*p+6*q+8));

E_F_BG = 0; k_h_BG = 0; k_If_BG = 0; k_Ir_BG = 0; % Placeholder for future extension

% Particle radius

x3 = [x_piv;x_piv;x_piv];
x5 = [x_piv;x_piv;x_piv;x_piv;x_piv];
R = sqrt(sum(y(p+q+1:6*p+6*q).*(162*x5+18))/(ro*pi*L*n));

% Total surface area

As = (n*2*pi*R*L)/1000;

if R-R0<=0
    R_ratio = 0;
else
    R_ratio = (R-R0)/R;
end

% L-HCM

rF_kin = zeros(2,1);    % Unregulated cellobiose uptake flux, mol/g-biom.s
rFE_kin = zeros(2,1);   % Unregulated inducive intracellular enzyme synthesis rates, 1/s
eF_rel = zeros(2,1);    % Relative intracellular enzyme levels
ROI = zeros(2,1);       % ROIs

for i = 1:2
    rF_kin(i) = kmax(i)*(y(2)/(K(i)+y(2)))*(1/(1+(y(1)/K_inh_G)))...
        *(1/(1+(y(6*p+6*q+12)/K_inh_etoh)));
    rFE_kin(i) = kE(i)*(y(2)/(K(i)+y(2)))*(1/(1+(y(1)/K_inh_G)))...
        *(1/(1+(y(6*p+6*q+12)/K_inh_etoh)));
    eF_rel(i) = y(i+6*p+6*q+10-1)/eF_max(i);
    ROI(i) = eF_rel(i)*rF_kin(i);
end

uF = zeros(2,1);        % Cybernetic variables
vF = zeros(2,1);        % Cybernetic variables
rF = zeros(2,1);        % Regulated cellobiose uptake flux, mol/g-biom.s

for i = 1:2
    uF(i) = ROI(i)/sum(ROI);
    vF(i) = ROI(i)/max(ROI);
    rF(i) = vF(i)*eF_rel(i)*rF_kin(i);
end

mu = Y2_B*rF(2);        % Specific biomass growth rate, 1/s

% Soluble polymers

dydt(1) = sum(n_exo(1,:)'.*k_h_exo.*y(p+q+1:2*p+2*q))...
    +sum(n_endo(1,7:end)'.*k_h_endo(7:end).*y(2*p+2*q+7:3*p+3*q))...
    +sum(n_endo(1,4:end)'.*k_h_endo(4:end).*y(3*p+3*q+4:4*p+4*q))...
    +sum(E_F_endo*n_endo(1,2:6)'.*k_hs_endo(2:6).*y(2:6))...
    +2*k_h_BG*E_F_BG*y(2)+k_Ir_exo(1)*y(6*p+6*q+5)+k_Ir_exo(1)*y(6*p+6*q+7)...
    +k_Ir_endo(1)*y(6*p+6*q+3)+k_Ir_BG*y(6*p+6*q+9)...
    -k_If_exo(1)*E_F_exo*y(1)-k_If_exo(1)*y(6*p+6*q+2)*y(1)...
    -k_If_endo(1)*E_F_endo*y(1)-k_If_BG*E_F_BG*y(1);

dydt(2) = sum(n_exo(2,:)'.*k_h_exo.*y(p+q+1:2*p+2*q))...
    +sum(n_endo(2,7:end)'.*k_h_endo(7:end).*y(2*p+2*q+7:3*p+3*q))...
    +sum(n_endo(2,4:end)'.*k_h_endo(4:end).*y(3*p+3*q+4:4*p+4*q))...
    +sum(E_F_endo*n_endo(2,3:6)'.*k_hs_endo(3:6).*y(3:6))...
    -k_hs_endo(2)*E_F_endo*y(2)-k_h_BG*E_F_BG*y(2)...
    +k_Ir_exo(2)*y(6*p+6*q+6)+k_Ir_exo(2)*y(6*p+6*q+8)+k_Ir_endo(2)*y(6*p+6*q+4)...
    -k_If_exo(2)*E_F_exo*y(2)-k_If_exo(2)*y(6*p+6*q+2)*y(2)-k_If_endo(2)*E_F_endo*y(2)...
    -sum(rF)*y(6*p+6*q+17);

for i = 3:5
    dydt(i) = sum(n_endo(i,7:end)'.*k_h_endo(7:end).*y(2*p+2*q+7:3*p+3*q))...
        +sum(n_endo(i,i+1:end)'.*k_h_endo(i+1:end).*y(3*p+3*q+i+1:4*p+4*q))...
        +sum(E_F_endo*n_endo(i,i+1:6)'.*k_hs_endo(i+1:6).*y(i+1:6))...
        -k_hs_endo(i)*E_F_endo*y(i);
end

for i = 6
    dydt(i) = sum(n_endo(i,i+1:end)'.*k_h_endo(i+1:end).*y(2*p+2*q+i+1:3*p+3*q))...
        +sum(n_endo(i,i+1:end)'.*k_h_endo(i+1:end).*y(3*p+3*q+i+1:4*p+4*q))...
        -k_hs_endo(i)*E_F_endo*y(i);
end

% Exo-enzyme-bound surface polymers

dydt(p+q+3) = sum(n_exo(3,3:end)'.*k_h_exo(3:end).*y(p+q+3:2*p+2*q));

for i = p+q+4:p+q+6
    dydt(i) = sum(n_exo(i-p-q,:)'.*k_h_exo.*y(p+q+1:2*p+2*q))...
        -(k_f_endo(i-p-q)*y(6*p+6*q+1)*y(i))...
        +(k_e_endo(i-p-q)*y(i-p-q+3*p+3*q));
end

for i = p+q+7:2*p+2*q-1
    dydt(i) = sum(n_exo(i-p-q,:)'.*k_h_exo.*y(p+q+1:2*p+2*q))...
        +(y(6*p+6*q+2)*k_f_exo(i-p-q)*y(i-p-q+4*p+4*q))...
        -(k_e_exo(i-p-q)*y(i))...
        -(k_f_endo(i-p-q)*y(6*p+6*q+1)*y(i))...
        +(k_e_endo(i-p-q)*y(i-p-q+3*p+3*q));

end

for i = 2*p+2*q
    dydt(i) = sum(n_exo(i-p-q,:)'.*k_h_exo.*y(p+q+1:2*p+2*q))...
        +(y(6*p+6*q+2)*k_f_exo(i-p-q)*y(i-p-q+4*p+4*q))...
        -(k_e_exo(i-p-q)*y(i))...
        -(k_f_endo(i-p-q)*y(6*p+6*q+1)*y(i))...
        +(k_e_endo(i-p-q)*y(i-p-q+3*p+3*q));
end

% Endo-enzyme-bound surface polymers

for i = 2*p+2*q+7:3*p+3*q
    dydt(i) = (k_f_endo(i-2*p-2*q)*y(6*p+6*q+1)*y(i-2*p-2*q+4*p+4*q))...
        -(k_e_endo(i-2*p-2*q)*y(i))...
        -(k_h_endo(i-2*p-2*q)*y(i));
end

% Exo-endo-bound surface polymers

for i = 3*p+3*q+4:4*p+4*q
    dydt(i) = (k_f_endo(i-3*p-3*q)*y(6*p+6*q+1)*y(i-3*p-3*q+p+q))...
        -(k_e_endo(i-3*p-3*q)*y(i))...
        -(k_h_endo(i-3*p-3*q)*y(i));
end

% Surface polymers

m = zeros(p+q,1);

for i = 4*p+4*q+7:5*p+5*q-1
    m(i-4*p-4*q) = (sum(n_endo(i-4*p-4*q,i-4*p-4*q+1:end)'...
        .*k_h_endo(i-4*p-4*q+1:end).*y(i-4*p-4*q+2*p+2*q+1:3*p+3*q))...
        +sum(n_endo(i-4*p-4*q,i-4*p-4*q+1:end)'...
        .*k_h_endo(i-4*p-4*q+1:end).*y(i-4*p-4*q+1+3*p+3*q:4*p+4*q))...
        -(k_f_exo(i-4*p-4*q)*y(6*p+6*q+2)*y(i))...
        +(k_e_exo(i-4*p-4*q)*y(i-4*p-4*q+p+q))...
        -(k_f_endo(i-4*p-4*q)*y(6*p+6*q+1)*y(i))...
        +(k_e_endo(i-4*p-4*q)*y(i-4*p-4*q+2*p+2*q)));
end

m(p+q) = (-(k_f_exo(p+q)*y(6*p+6*q+2)*y(5*p+5*q))...
    +(k_e_exo(p+q)*y(2*p+2*q))...
    -(k_f_endo(p+q)*y(6*p+6*q+1)*y(5*p+5*q))...
    +(k_e_endo(p+q)*y(3*p+3*q)));

for i = 4*p+4*q+7:5*p+5*q
    dydt(i) = m(i-4*p-4*q)-(sum(m.*(162*x_piv+18))...
        +sum(dydt(p+q+1:4*p+4*q).*(162*x3+18)))*R_ratio...
        *ppval(exp_ratio{i-4*p-4*q},R-R0);
end

% Internal inaccessible polymer

for i = 5*p+5*q+7:6*p+6*q
    dydt(i) = (sum(m.*(162*x_piv+18))...
        +sum(dydt(p+q+1:4*p+4*q).*(162*x3+18)))*R_ratio...
        *ppval(exp_ratio{i-5*p-5*q},R-R0);
end

% Surface-adsorbed endo-enzymes

dydt(6*p+6*q+1) = (k_ads_endo*E_F_endo*(As*tau_endo-y(6*p+6*q+1)))...
    -(k_des_endo*y(6*p+6*q+1))...
    -sum(y(6*p+6*q+1)*k_f_endo(7:end).*y(4*p+4*q+7:5*p+5*q))...
    +sum(k_e_endo(7:end).*y(2*p+2*q+7:3*p+3*q))...
    -sum(y(6*p+6*q+1)*k_f_endo(4:end).*y(p+q+4:2*p+2*q))...
    +sum(k_e_endo(4:end).*y(3*p+3*q+4:4*p+4*q));

% Surface-adsorbed exo-enzymes

dydt(6*p+6*q+2) = (k_ads_exo*E_F_exo*(As*tau_exo-y(6*p+6*q+2)))...
    -(k_des_exo*y(6*p+6*q+2))...
    -sum(y(6*p+6*q+2)*k_f_exo(7:end).*y(4*p+4*q+7:5*p+5*q))...
    +sum(k_e_exo(7:end).*y(p+q+7:2*p+2*q))...
    +k_Ir_exo(1)*y(6*p+6*q+7)+k_Ir_exo(2)*y(6*p+6*q+8)...
    -y(6*p+6*q+2)*sum(k_If_exo'.*y(1:2));

% Inhibited endo-enzymes

for i = 6*p+6*q+3:6*p+6*q+4
    dydt(i) = k_If_endo(i-6*p-6*q-3+1)*E_F_endo*y(i-6*p-6*q-3+1)...
        -k_Ir_endo(i-6*p-6*q-3+1)*y(i);
end

% Inhibited exo-enzymes

for i = 6*p+6*q+5:6*p+6*q+6
    dydt(i) = k_If_exo(i-6*p-6*q-5+1)*E_F_exo*y(i-6*p-6*q-5+1)...
        -k_Ir_exo(i-6*p-6*q-5+1)*y(i);
end

for i = 6*p+6*q+7:6*p+6*q+8
    dydt(i) = k_If_exo(i-6*p-6*q-7+1)*y(6*p+6*q+2)*y(i-6*p-6*q-7+1)...
        -k_Ir_exo(i-6*p-6*q-7+1)*y(i);
end

% Inhibited BG (placeholder for future extension)

for i = 6*p+6*q+9;
    dydt(i) = k_If_BG*E_F_BG*y(i-6*p-6*q-9+1)...
        -k_Ir_BG*y(i);
end

% Intracellular enzymes

for i = 6*p+6*q+10:6*p+6*q+11
    dydt(i) = aF(i-6*p-6*q-10+1)+uF(i-6*p-6*q-10+1)...
        *rFE_kin(i-6*p-6*q-10+1)-bF(i-6*p-6*q-10+1)*y(i)-mu*y(i);
end

% Metabolites

for i = 6*p+6*q+12:6*p+6*q+15
    dydt(i) = (Y1_met_lump(i-6*p-6*q-12+1)*rF(1)...
        +Y2_met_lump(i-6*p-6*q-12+1)*rF(2))*y(6*p+6*q+17);
end

% Cellulosome

dydt(6*p+6*q+16) = aE*y(6*p+6*q+17)+((Y1_E*rF(1)+Y2_E*rF(2))*y(6*p+6*q+17));

% Biomass

dydt(6*p+6*q+17) = mu*y(6*p+6*q+17);

% Output

ode = dydt;

end

%%-------------------------------------------------------------------------

function [C,CXS,CNS,CNXS,CFS,CS,CI,CT,CC,E_T_exo,E_T_endo,CB,etoh,lac,...
    form,ace,cellome,biom,cellulase_tot,cellu,cell_gluc_eq,eF_rel,R,...
    conv_rem,conv_cb,uF,vF,r_up,mu] = post_processing(t,y,p,q,MW_cellome,...
    mol_endo,mol_exo,mass_frac_cellulase,x_piv,eF_max,ro,L,n,mss,...
    cellobiose,kmax,K,K_inh_G,K_inh_etoh,Y2_B)

% Molar concentrations, mol/L

C = y(:,1:p+q);                     % Soluble products
CXS = y(:,p+q+1:2*p+2*q);           % CBH-bound complex
CNS = y(:,2*p+2*q+1:3*p+3*q);       % EG-bound complex
CNXS = y(:,3*p+3*q+1:4*p+4*q);      % CBH-EG-bound complex
CFS = y(:,4*p+4*q+1:5*p+5*q);       % Free un-bound surface polymers
CS = CXS+CNS+CNXS+CFS;              % Total surface polymers
CI = y(:,5*p+5*q+1:6*p+6*q);        % Internal polymers 
CT = C+CS+CI;                       % Total polymers
CC = CS+CI;                         % Insoluble polymers

E_T_endo = mol_endo*(y(:,6*p+6*q+16)/MW_cellome);   % Total endo-enzymes
E_T_exo = mol_exo*(y(:,6*p+6*q+16)/MW_cellome);     % Total exo-enzymes

E_S_endo = y(:,6*p+6*q+1);  % Surface-adsorbed endo-enzymes
E_S_exo = y(:,6*p+6*q+2);   % Surface-adsorbed exo-enzymes

CNI(:,1) = y(:,6*p+6*q+3);  % Inhibited endo-enzymes
CNI(:,2) = y(:,6*p+6*q+4);

CXI(:,1) = y(:,6*p+6*q+5)+y(:,6*p+6*q+7);   % Inhibited exo-enzymes
CXI(:,2) = y(:,6*p+6*q+6)+y(:,6*p+6*q+8);
  
E_F_endo = zeros(length(t),1);  % Free endo-enzymes
E_F_exo = zeros(length(t),1);   % Free exo-enzymes

for i = 1:length(t)
    E_F_endo(i) = E_T_endo(i)-E_S_endo(i)-sum(CNS(i,:))-sum(CNXS(i,:))-sum(CNI(i,:));
    E_F_exo(i) = E_T_exo(i)-E_S_exo(i)-sum(CXS(i,:))-sum(CXI(i,:))-sum(CNXS(i,:));
end

% Molecular weight of metabolites, g/mol

mw_met = [46.07;90.08;45.017;59.044];   

% Mass Concentrations, g/L

CB = y(:,2).*342;                       % Cellobiose
cellome = y(:,6*p+6*q+16);              % Cellulosome
biom = y(:,6*p+6*q+17);                 % Biomass
etoh = y(:,6*p+6*q+12).*mw_met(1);      % Ethanol
lac = y(:,6*p+6*q+13).*mw_met(2);       % Lactate
form = y(:,6*p+6*q+14).*mw_met(3);      % Formate
ace = y(:,6*p+6*q+15).*mw_met(4);       % Acetate

cellulase_tot = y(:,6*p+6*q+16).*mass_frac_cellulase;   % Total cellulase

cellu = zeros(length(t),1);     % Remaining insoluble cellulose

for i = 1:length(t)
    cellu(i) = sum(CC(i,:)'.*(162*x_piv+18));
end

cell_gluc_eq = cellu.*(180/162); % Remaining insoluble cellulose (g gluc/L)

% Relative intracellular enzyme levels

eF_rel(:,1) = y(:,6*p+6*q+10)/eF_max(1); 
eF_rel(:,2) = y(:,6*p+6*q+11)/eF_max(2);

% Transient particle radius, m

R = zeros(length(t),1);

for i = 1:length(t)
    R(i) = sqrt(sum(CC(i,:)'.*(162*x_piv+18))/(ro*pi*L*n));
end

% Conversion

conv_rem = zeros(length(t),1);  % Conversion of cellulose

for i = 1:length(t)
    conv_rem(i) = (mss-cellu(i))/mss;
end

conv_cb = zeros(length(t),1);  % Conversion of cellobiose 
for i = 1:length(t);
    conv_cb(i) = (cellobiose-CB(i))/cellobiose;
end

% Metabolic regulations

ROI = zeros(length(t),2);       % ROIs
rF_kin = zeros(length(t),2);    % Unregulated cellobiose uptake flux, mol/g-biom.s

for i = 1:length(t)
    rF_kin(i,1) = (kmax(1)*(y(i,2)/(K(1)+y(i,2)))...
        *(1/(1+(y(i,1)/K_inh_G)))*(1/(1+(y(i,6*p+6*q+12)/K_inh_etoh))));
    rF_kin(i,2) = (kmax(2)*(y(i,2)/(K(2)+y(i,2)))...
        *(1/(1+(y(i,1)/K_inh_G)))*(1/(1+(y(i,6*p+6*q+12)/K_inh_etoh))));
    ROI(i,1) = eF_rel(i,1)*rF_kin(i,1);
    ROI(i,2) = eF_rel(i,2)*rF_kin(i,2); 
end

uF = zeros(length(t),2);    % Cybernetic variables
vF = zeros(length(t),2);

for i = 1:length(t)
    uF(i,1) = ROI(i,1)/sum(ROI(i,:));
    uF(i,2) = ROI(i,2)/sum(ROI(i,:));
    vF(i,1) = ROI(i,1)/max(ROI(i,:));
    vF(i,2) = ROI(i,2)/max(ROI(i,:));
end

r_up = zeros(length(t),2);  % Regulated cellobiose uptake/fluxes, mol/g-biom.s

for i = 1:length(t)
    r_up(i,1) = vF(i,1)*eF_rel(i,1)*rF_kin(i,1);
    r_up(i,2) = vF(i,2)*eF_rel(i,2)*rF_kin(i,2);
end
    
mu = Y2_B.*r_up(:,2);       % Specific biomass growth rate, 1/s

end