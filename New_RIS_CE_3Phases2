clear all
clc
close all
%This code runs a montecarlo simulation for different noise powers
%TO DO:
%doesn't work for different power distribution between users
%----------------------------------------------------------------------------
%Definitions

M=16; %Number of BS antenna 
N_columns=8; %RIS columns
N_rows=8; %RIS rows
N=N_columns*N_rows; %number of RIS elements
K=10; %Number user antenna

c=300000000; %speed of light
f_c=2.4*10^9; %frequency
lambda=c/f_c; %wave length
d_BS=lambda/2; %BS inter antenna spacing
d_RIS=lambda/2; %RIS inter element spacing

d_RIS=50; %vertical distance between BS and RIS / RIS vertical position
h_BS=10; %BS height
h_RIS=20; %RIS height
d_1=sqrt(d_RIS^2+(h_BS-h_RIS)^2); %distance between BS and RIS
d_user=40*ones(1,K); %user vertical position relative to BS
d_k=sqrt(d_user.^2+h_BS^2);%distance between BS and user k
d_2k=sqrt((d_user-d_RIS).^2+h_RIS^2);%distance between RIS and user k
theta_LOS=atan(abs(h_BS-h_RIS)/d_RIS); %LOS elevation angle

p_k=ones(K,1);%transmission power for indexed user (Kx1) (unused)
P=diag(p_k,0); %(unused)
P_C=1; %this is the only used power, same for all users. The code can be
%changed to work with other power allocations using p_k

tau_1=0.1; %coherence period
ttau=10*10^(-6); %symbol period


%channel
R_RISk=zeros(N,N,K);%RIS correlation matrix (NxNxK)
R_BSk=zeros(M,M,K);%BS correlation matrix (MxMxK)
for k=1:K
    R_RISk(:,:,k)=eye(N);%RIS correlation matrix (NxNxK)
    R_BSk(:,:,k)=eye(M);%BS correlation matrix (MxMxK)
end
alpha_1=2.2;
alpha_2=3.67;
C_1=25;
C_2=30;
beta_2k=10^((-C_2)/10)./(d_2k.^alpha_2);%pathloss from RIS to user k (Kx1)
beta_dk=10^((-C_2)/10)./(d_k.^alpha_2);%pathloss from BS to user k (Kx1)
beta_1=10^((-C_1)/10)/(d_1^alpha_1); %H_1 pathloss

%iteration variables, for noise and monte carlo repitions
k_b=1.38*10^(-23); %boltzman constant
Tem_k=300; %temperature in Kelvin
band=1/ttau; %largura de banda
sigma_sq=k_b*Tem_k*band; %
sigma2db=pow2db(sigma_sq); %noise squared in db
sigma=sqrt(sigma_sq);
%sigma2db=-180;MC.mcPLot

noise_range=-175:5:-80;%'Noise'
distance_range=0:10:100; % 'Distance'
BS_range=18:20:128; %'BS'
RIS_range=4:3:16; % 'RIS'

channel = Channel(M,N_columns,N_rows,K,d_user,sigma2db);
channel = channel.setPropagation_loss(40);
MC=MonteCarlo(100,channel); %doesnt work for 1
MC=MC.flagPlotAvgCH3P(0);
MC=MC.flagPlotp3pCH3P(0);
MC=MC.flagNoiseSNR(0);
MC=MC.flagLSestimation(1);
MC=MC.flagLambdaEstimation(0);
MC=MC.fixSubphases(0);
MC=MC.fixNsymbols(0);
MC=MC.fixP1symbols(0);
MC=MC.fixP2symbols(0);
MC=MC.fixP3symbols(0); %For this to work properly, M=N and T_3=a*(K-1) for some natural a
MC=mcProcess(MC,noise_range,'all','Noise');
MC.mcPlotCapacity();
MC.mcPlot2();
%MC.plot_pathloss();
