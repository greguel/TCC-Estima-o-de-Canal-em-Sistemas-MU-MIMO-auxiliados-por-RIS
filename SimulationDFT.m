classdef SimulationDFT
        properties
        channel

        g
        channel_RIS_Users
        channel_BS_RIS
        channel_BS_Users

        phi
        effective_channel
        
        S
        s %symbols in phase 1
        T_s %number of symbols in phase 1
        tau %coherence period of phase 1
        h_dkest %direct channel estimation
        g_est

        h_dkls
        g_ls

        stats_DFT1
        stats_DFT2
        stats_DFTLS1
        stats_DFTLS2

        psi_dk
        psi_nk

        user_rate_DFT
        user_rate_LS
        user_rate_PE
    end
   
    
    methods
        %constructor
        function obj=SimulationDFT(channel)
            obj.channel=channel;
            obj.g=[];
            obj.effective_channel=0;
            obj.phi=[];
            obj.S=channel.RIS_elements+1;
            obj.T_s=channel.Users;
            obj.tau=0;
            obj.s=[];
            obj.h_dkest=[];
            obj.channel_RIS_Users=[];
            obj.channel_BS_RIS=[];
            obj.channel_BS_Users=[];
            obj.h_dkls=[];
            obj.g_ls=[];


            K=channel.Users;

            obj.stats_DFT1=Statistics();
            obj.stats_DFT2=repmat(Statistics, 1, K);
            obj.stats_DFTLS1=Statistics();
            obj.stats_DFTLS2=repmat(Statistics, 1, K);
            obj.psi_dk=[];
            obj.psi_nk=[];
            obj.user_rate_DFT=0;
            obj.user_rate_LS=0;
            obj.user_rate_PE=0;
        end


        function obj=channelEstimation(obj) %channel estimation
            channel=obj.channel;
            N=channel.RIS_elements;
            K=channel.Users;
            M=channel.BS_antennas;
            S=obj.S;
            ttau=channel.symbol_period;
            P_C=channel.power_Users;
            sigma=channel.noise_abs;
            beta_dk=channel.pathloss_BS_Users;
            beta_2k=channel.pathloss_RIS_Users;
            R_RISk=channel.RIS_correlation_matrix;
            R_BSk=channel.BS_correlation_matrix;
            H_1=obj.channel_BS_RIS;

            %Definitions, space allocation or just unused variables in code
            alpha=ones(N,1); %vector of RIS element gains (Nx1)
            theta=ones(N,1); %vector of RIS element phase shifts (Nx1)
            Theta=diag(alpha.*exp(1i*theta)); %phase shift matrix (NxN)

            v=diag(Theta); %equivalent reflective beamforming vector unused

            alpha_sn=ones(N,S); %reflect beamform gain at element n at subphase s unused
            theta_sn=ones(N,S); %reflect beamform gain at element n at subphase s unused
            v_s=alpha_sn.*exp(1i*theta_sn); %reflect beamform applied during subphase s (NxS)

            for ii=1:S%V DFT matrix for estimation
                for k=1:N+1
                    V(ii,k)=exp(-1i*2*pi*(ii-1)*(k-1)/S);
                end
            end
            Vb=kron(V,eye(M)); %size (MSxM(N+1))
            for k=1:S
                v_s(:,k)=transpose(V(k,2:end)); %this is the used subphase RIS phase matrix
            end
            %--------------------------------------------
            %signal model
            
            obj=obj.make_orthogonal();
            s=obj.s;
            T_s=obj.T_s;

            tau_s=T_s*ttau; %tau_c/S subphase traning period
            tau_c=tau_s*S; %uplink traning period
            %tau_d=tau-tau_c; %downlink transmission period

            s=s;%symbol for indexed user, s_k is s_k(1:k) (T_SxK)

            x=s*sqrt(P_C); %needs to be scaled by power

            for ii=1:S
                obj.phi=v_s(:,ii);
                [obj,Y_s(:,:,ii)]=send_signal(obj,T_s,s,1,P_C,sigma);
            end %received (MxT_sxS) signal

            for ii=1:S
                for k=1:K
                    r_sk(:,ii,k)=Y_s(:,:,ii)*s(k,:)'/(sqrt(P_C)*T_s);
                end
            end %(MxSxK)signal correlated by x, contains channel information and noise

            r_k=reshape(r_sk(:,1,:),M,K);
            for ii=2:S
                r_k= cat(1,r_k,reshape(r_sk(:,ii,:),M,K));
            end %the correllated signal concatenated, that is, r_sk is stacked to become (MSxK)
            Vb_pinv = kron(pinv(V),eye(M));
            for k=1:K
                rt_k(:,k)=Vb_pinv*r_k(:,k); %Last changed
            end %final matrix before estimation, found by using the pseudo inverse of Vb to "remove it"
            %all the channels are stacked here and have noise added, M direct channels
            % and MN channels that go through the RIS(M+MNxK). This is also
            % the LS estimator
            obj.h_dkls=rt_k(1:M,:);
            for k=1:K
                for ii=2:N+1
                    obj.g_ls(:,ii-1,k)=rt_k(M*(ii-1)+1:M*ii,k);
                end
            end

            %deconcatenating rt_k, this is the rt_1k, rt_2k, etc. matrix
            for k=1:K
                for ii=1:N+1
                    rt_ik(:,ii,k)=rt_k(M*(ii-1)+1:M*ii,k);
                end
            end %this is the destacked and organized rt_k, the first column is the
            % the direct channels, the other columns the RIS channels (MxN+1xK)

            %------------------------------------------------------------------------------
            %Estimation
            for k=1:K %just used for calculations
                Q_dk(:,:,k)=inv(beta_dk(k)*R_BSk(:,:,k)+sigma^2/(S*P_C*T_s)*eye(M));
            end

            %estimate of h_dk
            for k=1:K
                obj.h_dkest(:,k)=beta_dk(k)*R_BSk(:,:,k)*Q_dk(:,:,k)*rt_ik(:,1,k);
            end

            for k=1:K %just used for calculations
                for ii=1:N
                    Q_nk(:,:,ii,k)=inv(R_RISk(ii,ii,k)*beta_2k(k)*H_1(:,ii)*H_1(:,ii)'+sigma^2/(S*P_C*T_s)*eye(M));
                end
            end
            %estimate of H_0k
            for k=1:K
                for ii=1:N
                    obj.g_est(:,ii,k)=R_RISk(ii,ii,k)*beta_2k(k)*H_1(:,ii)*H_1(:,ii)'*Q_nk(:,:,ii,k)*rt_ik(:,ii+1,k);
                end
            end

            %-----------------------------------------------------------------------------
            %analytical error variances
            for k=1:K
                obj.psi_dk(:,:,k)=beta_dk(k)*R_BSk(:,:,k)-beta_dk(k)^2*R_BSk(:,:,k)*Q_dk(:,:,k)*R_BSk(:,:,k);
            end

            for k=1:K
                for ii=1:N
                    obj.psi_nk(:,:,ii,k)=beta_2k(k)*R_RISk(ii,ii,k)*H_1(:,ii)*H_1(:,ii)'-R_RISk(ii,ii,k)*R_RISk(ii,ii,k)'*beta_2k(k)^2*H_1(:,ii)*H_1(:,ii)'*Q_nk(:,:,ii,k)*H_1(:,ii)*H_1(:,ii)';
                end
            end
            obj = obj.channelCapacity;
        end
        

        function obj=simulationStats(obj) %calculates the statistics for the channel estimation protcol
            channel=obj.channel;
            h_dkest=obj.h_dkest;
            h_dkls=obj.h_dkls;
            g=obj.g;
            g_est=obj.g_est;
            g_ls=obj.g_ls;
            H_1=obj.channel_BS_RIS;
            h_dk=obj.channel_BS_Users;
            M=channel.BS_antennas;
            K=channel.Users;
            N=channel.RIS_elements;
            beta_1=channel.pathloss_BS_RIS;
            beta_2k=channel.pathloss_RIS_Users;
            beta_dk=channel.pathloss_BS_Users;
            sigma=channel.noise_abs;
            S=obj.S;
            P_C=channel.power_Users;
            T_s=obj.T_s;
            ttau=channel.symbol_period;

            obj.stats_DFT1=stats(obj.stats_DFT1,h_dk,h_dkest,M*K);
            obj.stats_DFTLS1=stats(obj.stats_DFTLS1,h_dk,h_dkls,M*K);
            
            for k=1:K
                obj.stats_DFT2(k)=stats(obj.stats_DFT2(k),g(:,:,k),g_est(:,:,k),M*N);
                obj.stats_DFTLS2(k)=stats(obj.stats_DFTLS2(k),g(:,:,k),g_ls(:,:,k),M*N);
            end

            %theoretical NMSE
            temp2=zeros(1,K);
            temp=zeros(N,K);
            for ii=1:K %averaging for beta_dk
                temp2(1,ii)=(sigma^2/(S*P_C*T_s))/(beta_dk(ii)+sigma^2/(S*P_C*T_s));
                for jj=1:N %used for finding the theoretical NMSE
                    A=beta_2k(1)*H_1(:,jj)*H_1(:,jj)'+sigma^2/(S*P_C*T_s)*eye(M);
                    precond_M = diag(diag(A));
                    M_inv = pinv(precond_M);
                    A_prime = M_inv * A;


                    b_prime = M_inv * eye(size(A,1));
                    A_inv = zeros(size(A,1));

                    for i = 1:size(A,1)
                        A_inv(:, i) = A_prime \ b_prime(:, i); % Solve A' * x = b'
                    end

                    temp(jj,ii)=1/(M*beta_1*beta_2k(1))*(beta_2k(1)*trace(H_1(:,jj)*H_1(:,jj)')-beta_2k(1)^2*trace(H_1(:,jj)*H_1(:,jj)'*A_inv*H_1(:,jj)*H_1(:,jj)'));
                end
            end
            obj.stats_DFT1.t_nmse=mean(temp2);
            %t_nmse_nk(k)=(sigma(k)^2/(S*P_C*T_s))/(M*beta_2k(1)+sigma(k)^2/(S*P_C*T_s));
            %%this is for when H_1 has rank 1

            obj.stats_DFT2(1,1).t_nmse=mean(abs(temp),'all'); %NMSE for generic H_1
        end

        
        function obj=find_effective_channel(obj) %finds the effective channel, used before sending a signal
            channel=obj.channel;
            K=channel.Users;
            h_dk=obj.channel_BS_Users;
            H_1=obj.channel_BS_RIS;
            h_2k=obj.channel_RIS_Users;
            phi=obj.phi;

            for k=1:K
                H(:,k)=h_dk(:,k)+H_1*diag(h_2k(:,k))*phi; %effective channel
            end
            obj.effective_channel=H;
        end

        
        function obj=make_orthogonal(obj) %makes orthogonal singals
            T_s=obj.T_s;
            K=obj.channel.Users;
            obj.s=hadamard(2^ceil(log2(T_s))); %generates orthogonal pilot signals
            for i=1:2^ceil(log2(T_s))-K %removes unecessary rows
                obj.s(2^ceil(log2(T_s))-i,:)=[];
            end
            obj.T_s=2^ceil(log2(T_s)); %number of training symbols in first phase
        end

        
        function [obj,y]=send_signal(obj,T,s,ttau,P_C,sigma) %sends a signal through RIS channel
            obj=find_effective_channel(obj);
            M=obj.channel.BS_antennas;
            H=obj.effective_channel;
            x_ki=s*sqrt(P_C); %signal for user k at time i

            n=(randn(M,T) + 1i*randn(M,T))*sigma/sqrt(2);
            %n=(randn(M,I) + 1i*randn(M,I))*sigma/sqrt(2); %phase noise
            y=H*x_ki+n; %signal model
        end
        

        function obj=channel_avg(obj) %calculates the channel averages needed for CE
            for j=1:500
                channel=obj.channel;
                obj=obj.makeChannels();
                g=obj.g;
                lamb=obj.lamb;
                K=channel.Users;
                N=channel.RIS_elements;
                M=channel.BS_antennas;
                T_3=(K-1)*ceil(N/M);
                

                temp1(:,:,j)=g(:,:,1)'*g(:,:,1);
                for ii=1:T_3
                    phase_3(ii).k=ceil(ii/ceil(N/M))+1;
                    phase_3(ii).first=(ii-floor(ii/ceil(N/M))*ceil(N/M)-1)*M+1:(ii-floor(ii/ceil(N/M))*ceil(N/M)-1)*M+M;
                    phase_3(ii).conditional=(floor(ii/ceil(N/M))~=ii/ceil(N/M));
                    phase_3(ii).second=(ceil(N/M)-1)*M+1:N;
                    if phase_3(ii).conditional==0
                        phase_3(ii).delta=phase_3(ii).second;
                    else
                        phase_3(ii).delta= phase_3(ii).first;
                    end
                    for jj=1:size(phase_3(ii).delta,2)
                        phase_3(ii).Lamb(jj)=lamb(phase_3(ii).delta(jj),phase_3(ii).k);
                    end
                    temp2(ii).array(:,:,j)=phase_3(ii).Lamb'*phase_3(ii).Lamb;
                end
            end
            C_BI=mean(temp1,3);
            %conditioning
            A=max(C_BI,[],'all');
            cond_constant=10^(log10(abs(max(A(:))))*0.003);
            cond_matrix=cond_constant*eye(size(C_BI));
            obj.C_BI=C_BI+cond_constant;

            for ii =1:T_3
                %conditioning
                C_lamb{ii}=mean(temp2(ii).array,3);
                A=max(C_lamb{ii},[],'all');
                cond_constant=10^(log10(abs(max(A(:)))*0.3));
                cond_matrix=cond_constant*eye(size(C_lamb{ii}));
                obj.C_lamb{ii}=C_lamb{ii}+cond_matrix;
            end
            obj=obj.makeChannels;
        end


        function obj=updateChannel(obj,channel) %creates a new channel instance
            obj.channel=channel;
        end


        function obj = makeChannels(obj)
            %LoS channel between the BS and RIS. This assumed to be known
            % H_1(:,n) is the channel between the BS and element n (MxN),
            % should be rank 1, but the code should work with other channels
            channel=obj.channel;
            AoD_BS=[0,channel.elevation_BS_RIS]; %BS AOD azimuth/elevation
            AoA_RIS=[0,-channel.elevation_BS_RIS]; %RIS AOA, azimuth/elevation
            M=channel.BS_antennas;
            N=channel.RIS_elements;
            K=channel.Users;
            d_BS=channel.distance_BS_elements;
            d_RIS=channel.distance_RIS_elements;
            lambda=channel.wave_length;
            N_rows=channel.RIS_rows;
            N_columns=channel.RIS_columns;
            beta_1=channel.pathloss_BS_RIS;
            beta_2k=channel.pathloss_RIS_Users;
            beta_dk=channel.pathloss_BS_Users;
            R_BSk=channel.BS_correlation_matrix;
            R_RISk=channel.RIS_correlation_matrix;
            prop_l = channel.BS_Users_propogation_loss;
            rice_factor = channel.rice_factor;

            %array response
            a_H1=exp(-1i*2*pi*sin(AoD_BS(1))*cos(AoD_BS(2)).*(0:M-1)*d_BS/lambda);%BS array response ULA
            b_H1=kron(exp(-1i*2*pi*sin(AoA_RIS(2)).*(0:N_rows-1)*d_RIS/lambda),...
                exp(-1i*2*pi*sin(AoA_RIS(1))*cos(AoA_RIS(2)).*(0:N_columns-1)*d_RIS/lambda));%RIS array response UPA
            H_1=sqrt(beta_1)*transpose(a_H1)*b_H1;%scales by pathloss



            z_k=(randn(N,K) + 1i*randn(N,K))/sqrt(2);%fast fadding vector from RIS to user k (NxK)
            z_dk=(randn(M,K) + 1i*randn(M,K))/sqrt(2);%fast fadding vector from BS to user k (MxK)

            h_2k = zeros(N, K);
            h_dk = zeros(M, K);
            for jj=1:K
                h_2k(:,jj)=sqrt(beta_2k(jj))*sqrtm(R_RISk(:,:,jj))*z_k(:,jj); %RIS-UEk channel link
                % Nx1 for each user k (NxK), for user k h_2k(:,k)
                h_dk(:,jj)=sqrt(beta_dk(jj))*sqrtm(R_BSk(:,:,jj))*z_dk(:,jj); %BS-UE direct channel,
                % Mx1 for each user k (MxK),for user k h_dk(:,k)
            end

            H_1=H_1+sqrt(beta_1/db2pow(rice_factor))*(randn(size(H_1)) + 1i*randn(size(H_1)))/sqrt(2); %Rice factor

            obj.channel_BS_RIS=H_1; %Channel 1

            %h_2k=h_2k+10^(-6)*randn(size(h_2k));

            obj.channel_RIS_Users=h_2k;
            obj.channel_BS_Users=h_dk;

            for k=1:K
                obj.g(:,:,k)=H_1*diag(h_2k(:,k));
            end
            K=obj.channel.Users;
            N=obj.channel.RIS_elements;
            g=obj.g();
        end

        function obj=setSubphases(obj,S)
            obj.S=S;
        end

        function obj = setNsymbols(obj,T_s)
            obj.T_s=T_s;
        end

        function obj = channelCapacity(obj)
            H_0_est = obj.g_est(:,:,1);
            h_d_est = obj.h_dkest(:,1);
            channel = obj.channel;
            P = channel.power_Users;
            sigma = channel.noise_abs;
            h_d_ls = obj.h_dkls(:,1);
            H_0_ls = obj.g_ls(:,:,1);
            H_0 = obj.g(:,:,1);
            h_d = obj.channel_BS_Users(:,1);
            

            v = exp(1j * phase(H_0_est'*h_d_est));
            obj.phi = v;
            obj = obj.find_effective_channel;
            H = obj.effective_channel(:,1);

            H_est = h_d_est+H_0_est*v;
            w_k = H_est/norm(H_est);

            SNR = P*abs(w_k'*(H))^2/sigma^2;
            user_rate_DFT = log2(1+SNR)*(1-obj.T_s*obj.S*channel.symbol_period/channel.coherence_period)/channel.symbol_period;

            v = exp(1j * phase(H_0_ls'*h_d_ls));
            obj.phi = v;
            obj = obj.find_effective_channel;
            H = obj.effective_channel(:,1);

            H_ls = h_d_ls+H_0_ls*v;
            w_k = H_ls/norm(H_ls);

            SNR = P*abs(w_k'*(H))^2/sigma^2;
            user_rate_LS = log2(1+SNR)*(1-obj.T_s*obj.S*channel.symbol_period/channel.coherence_period)/channel.symbol_period;

            v = exp(1j * phase(H_0'*h_d));
            obj.phi = v;
            obj = obj.find_effective_channel;
            H = obj.effective_channel(:,1);

            w_k = H/norm(H);

            SNR = P*abs(w_k'*(H))^2/sigma^2;
            user_rate_PE = log2(1+SNR)*(1-obj.S*obj.T_s*channel.symbol_period/channel.coherence_period)/channel.symbol_period;

            obj.user_rate_PE = user_rate_PE;
            obj.user_rate_DFT = user_rate_DFT;
            obj.user_rate_LS = user_rate_LS;
            %rate = (1-estimation_period/coherence_period)*capacity;
        end
    end
end

