classdef Simulation3P
    properties
        channel

        g
        lamb
        channel_RIS_Users
        channel_BS_RIS
        channel_BS_Users

        phi
        effective_channel

        s_1 %symbols in phase 1
        T_1 %number of symbols in phase 1
        tau_1 %coherence period of phase 1
        h_dkest %direct channel estimation
        p1_var
        p1_avg

        s_2 %symbols in phase 1
        T_2 %number of symbols in phase 1
        tau_2 %coherence period of phase 1
        g_est %direct channel estimation
        g_estp %direct channel estimation with perfect estimation in phase 2
        psi_2
        Phi_2
        p2_var
        p2_avg

        s_3 %symbols in phase 1
        T_3 %number of symbols in phase 1
        tau_3 %coherence period of phase 1
        lamb_est %direct channel estimation
        lamb_estp
        Gp
        p3_var
        p3_avg
        psi_3

        stats_p1
        stats_p2
        stats_p3
        stats_p3p
        stats_p3l
        stats_p3lp

        C_BI
        C_lamb

        user_rate_3P
        user_rate_DC
    end
   
    
    methods
        %constructor
        function obj=Simulation3P(channel)
            obj.channel=channel;
            obj.g=[];
            obj.lamb;
            obj.effective_channel=0;
            obj.phi=[];
            obj.T_1=channel.Users;
            obj.tau_1=0;
            obj.s_1=[];
            obj.h_dkest=[];
            obj.s_2=[];
            obj.T_2=channel.RIS_elements; %second phase number of symbols 
            obj.tau_2=0; 
            obj.g_est=[];
            obj.g_estp=[];
            obj.s_3=[];
            
            obj.tau_3=0;
            obj.lamb_est=[]; 
            obj.lamb_estp=[];
            obj.lamb=[];
            obj.C_BI=0;
            obj.C_lamb=[];
            obj.channel_RIS_Users=[];
            obj.channel_BS_RIS=[];
            obj.channel_BS_Users=[];
            obj.psi_2=[];
            obj.Phi_2=[];
            obj.psi_3=[];
            obj.Gp=[];

            K=channel.Users;
            N=channel.RIS_elements;
            M=channel.BS_antennas;
            
            obj.T_3=(K-1)*ceil(N/M);
            obj.stats_p1=Statistics();
            obj.stats_p2=Statistics();
            obj.stats_p3=repmat(Statistics, 1, K-1);
            obj.stats_p3p=repmat(Statistics, 1, K-1);
            obj.stats_p3l=Statistics();
            obj.stats_p3lp=Statistics();

            obj.user_rate_3P=0;
            obj.user_rate_DC=0;
        end


        function obj=firstPhase(obj) %first phase
        channel=obj.channel;
        N=channel.RIS_elements;
        K=channel.Users;
        M=channel.BS_antennas;
        ttau=channel.symbol_period;
        P_C=channel.power_Users;
        sigma=channel.noise_abs;
        beta_dk=channel.pathloss_BS_Users;

        obj.phi=zeros(N,1); %RIS off

        obj=make_orthogonal(obj);
        
        T_1=obj.T_1;
        s_1=obj.s_1;%symbol for indexed user, s_k is s_k(1:k) (T_SxK)
        
        tau_1=T_1*ttau;
        obj.tau_1=tau_1; %phase traning period
        
        [obj,y_1i]=send_signal(obj,T_1,s_1,1,P_C,sigma);

        obj.h_dkest=beta_dk(1)*sqrt(P_C)*y_1i*s_1'/(beta_dk(1)*P_C*T_1+sigma^2); %direct channel estimation

        end


        function obj=secondPhase(obj)
            channel=obj.channel;
            N=channel.RIS_elements;
            K=channel.Users;
            M=channel.BS_antennas;
            ttau=channel.symbol_period;
            P_C=channel.power_Users;
            sigma=channel.noise_abs;
            beta_1=channel.pathloss_BS_RIS;
            T_1=obj.T_1;

            
            T_2=obj.T_2;

            tau_2=T_2*ttau; %phase length
            Phi_2=zeros(N,T_2);
            for k=1:T_2% DFT matrix for phi
                for ii=1:N
                    Phi_2(ii,k)=exp(-1i*2*pi*(ii-1)*(k-1)/T_2);
                end
            end
            obj.Phi_2=Phi_2;

            s_2=zeros(K-1,T_2);
            s_2=[ones(1,T_2);s_2];
            obj.s_2=s_2;
            x_2ki=s_2*sqrt(1*P_C);%symbol is 1 for user 1, zero for other users
            n=(randn(M,T_2) + 1i*randn(M,T_2))*sigma/sqrt(2);
            for ii=1:T_2
                obj.phi=Phi_2(:,ii);
                [obj,y_2i(:,ii)]=send_signal(obj,1,s_2(:,ii),1,P_C,sigma);
            end

            y_2ib=y_2i-obj.h_dkest*x_2ki;
            psi_2=P_C*M*beta_1*sigma^2*s_2(1,:)'*s_2(1,:)/(beta_1*P_C*T_1+sigma^2)+M*sigma^2*eye(T_2);
            C_BI=obj.C_BI; 
            obj.psi_2=psi_2;

            %g_est(:,:,1)=1/(sqrt(ttau*P_C)*T_2)*y_2ib*Phi_2'; %g estimation for user 1 no noise
            obj.g_est(:,:,1)=sqrt(P_C)*y_2ib*inv(psi_2)*Phi_2'*pinv(P_C*Phi_2*inv(psi_2)*Phi_2'+pinv(C_BI));
        end


        function obj=thirdPhase(obj)
            channel=obj.channel;
            N=channel.RIS_elements;
            K=channel.Users;
            M=channel.BS_antennas;
            ttau=channel.symbol_period;
            P_C=channel.power_Users;
            sigma=channel.noise_abs;
            beta_dk=channel.pathloss_BS_Users;
            T_1=obj.T_1;
            R_RISk=channel.RIS_correlation_matrix;
            R_BSk=channel.BS_correlation_matrix;
            g_est=obj.g_est;
            g=obj.g;
            h_dkest=obj.h_dkest;

            T_3=obj.T_3;

            n=(randn(M,T_3) + 1i*randn(M,T_3))*sigma/sqrt(2);
            psi_3=beta_dk(1)*P_C*sigma^4/(beta_dk(1)*P_C*T_1+sigma^2)*R_BSk(:,:,1)+(beta_dk(1)*P_C)^2*T_1*sigma^2/(beta_dk(1)*P_C*T_1+sigma^2)*eye(M)+sigma^2*eye(M);
            obj.psi_3=psi_3;

            for ii=1:(K-1)*ceil(N/M)
                phase_3(ii).k=ceil(ii/ceil(N/M))+1;
                phase_3(ii).first=(ii-floor(ii/ceil(N/M))*(ceil(N/M))-1)*M+1:(ii-floor(ii/ceil(N/M))*ceil(N/M)-1)*M+M;
                phase_3(ii).conditional=(floor(ii/ceil(N/M))~=ii/ceil(N/M));
                phase_3(ii).second=(ceil(N/M)-1)*M+1:N;
                if phase_3(ii).conditional==0
                    phase_3(ii).delta=phase_3(ii).second;
                else
                    phase_3(ii).delta= phase_3(ii).first;
                end
                for jj=1:size(phase_3(ii).delta,2)
                    phase_3(ii).G(:,jj)=g_est(:,phase_3(ii).delta(jj),1);
                    phase_3(ii).Gp(:,jj)=g(:,phase_3(ii).delta(jj),1);
                    phase_3(ii).Lamb(jj)=obj.lamb(phase_3(ii).delta(jj),phase_3(ii).k);
                end
                obj.Gp=phase_3(ii).Gp;
                phase_3(ii).Clamb=obj.C_lamb{ii};
            end
            

            for ii=1:T_3
                if T_3~=(K-1)*ceil(N/M)
                    Symbolflag=ceil(ii/(T_3/(K-1)));
                else
                    Symbolflag=ii;
                end
                for k=2:K
                    if sum(ismember(phase_3(Symbolflag).k,k))>=1
                        s_3(k,ii)=1;
                    else
                        s_3(k,ii)=0;
                    end
                end
                for jj=1:N
                    if sum(ismember(phase_3(Symbolflag).delta,jj))>=1
                        Phi_3(ii,jj)=1;
                    else
                        Phi_3(ii,jj)=0;
                    end
                end
                obj.phi=Phi_3(ii,:)';
                x_3ki(:,ii)=s_3(:,ii)*sqrt(P_C); 

                [obj,y_3i(:,ii)]=send_signal(obj,1,s_3(:,ii),1,P_C,sigma);
                
                %should be equal to
                %y_3ib(:,ii)=sqrt(P_C)*phase_3(ii).G*transpose(phase_3(ii).Lamb)
                y_3ib(:,ii)=y_3i(:,ii)-h_dkest*x_3ki(:,ii);
                %conditioning:
            end

            y_3ibavg=[];
            ii=1;
            if T_3~=(K-1)*ceil(N/M)
                while ii<=T_3
                    avg_columns = mean(y_3ib(:, ii:ii+T_3/(K-1)-1), 2);
                    y_3ibavg = [y_3ibavg avg_columns];

                    ii = ii + T_3/(K-1);
                end
            else
                y_3ibavg=y_3ib;
            end

            for ii=1:(K-1)*(ceil(N/M))
                %preconditioning
                A=P_C*phase_3(ii).G'*pinv(psi_3)*phase_3(ii).G+pinv(phase_3(ii).Clamb);
                precond_M = diag(diag(A));
                M_inv = inv(precond_M);
                A_prime = M_inv * A; 

              
                b_prime = M_inv * eye(size(A,1)); 
                A_inv = zeros(size(A,1)); 

                for i = 1:size(A,1)
                    A_inv(:, i) = A_prime \ b_prime(:, i); % Solve A' * x = b'
                end
                
                B=P_C*phase_3(ii).Gp'*pinv(psi_3)*phase_3(ii).Gp+pinv(phase_3(ii).Clamb);
                precond_M = diag(diag(B));
                M_inv = inv(precond_M);
                B_prime = M_inv * B;


                b_prime = M_inv * eye(size(B,1));
                B_inv = zeros(size(B,1));

                for i = 1:size(B,1)
                    B_inv(:, i) = B_prime \ b_prime(:, i); % Solve A' * x = b'
                end

               

                phase_3(ii).Lambest=sqrt(P_C)*A_inv*phase_3(ii).G'*inv(psi_3)*y_3ibavg(:,ii);
                phase_3(ii).Lambestp=sqrt(P_C)*B_inv*phase_3(ii).Gp'*inv(psi_3)*y_3ibavg(:,ii);
                for jj=1:size(phase_3(ii).Lambest,1)
                    obj.lamb_est(phase_3(ii).delta(jj),phase_3(ii).k)=phase_3(ii).Lambest(jj);
                    obj.lamb_estp(phase_3(ii).delta(jj),phase_3(ii).k)=phase_3(ii).Lambestp(jj);
                end
            end
            for k=2:K
                for ii=1:N
                    obj.g_est(:,ii,k)=obj.lamb_est(ii,k)*g_est(:,ii,1);
                    obj.g_estp(:,ii,k)=obj.lamb_estp(ii,k)*g(:,ii,1);
                end
            end

            obj = obj.channelCapacity;
        end
        
        
        function obj=simulationStats(obj) %calculates the statistics for the channel estimation protcol
            channel=obj.channel;
            h_dkest=obj.h_dkest;
            g=obj.g;
            g_est=obj.g_est;
            g_estp=obj.g_estp;
            lamb=obj.lamb;
            lamb_est=obj.lamb_est;
            lamb_estp=obj.lamb_estp;
            h_dk=obj.channel_BS_Users;
            M=channel.BS_antennas;
            K=channel.Users;
            N=channel.RIS_elements;
            beta_dk=channel.pathloss_BS_Users;
            P_C=channel.power_Users;
            T_1=obj.T_1;
            sigma=channel.noise_abs;
            ttau=channel.symbol_period;
            Phi_2=obj.Phi_2;
            psi_2=obj.psi_2;
            psi_3=obj.psi_3;
            C_BI=obj.C_BI;

            obj.stats_p1=stats(obj.stats_p1,h_dk,h_dkest,M*K);
            obj.stats_p2=stats(obj.stats_p2,g(:,:,1),g_est(:,:,1),M*N);
            obj.stats_p3l=stats(obj.stats_p2,lamb(:,2:end),lamb_est(:,2:end),N*(K-1));
            obj.stats_p3lp=stats(obj.stats_p2,lamb(:,2:end),lamb_estp(:,2:end),N*(K-1));
            for k=2:K
                obj.stats_p3(k-1)=stats(obj.stats_p3(k-1),g(:,:,k),g_est(:,:,k),M*N);
                obj.stats_p3p(k-1)=stats(obj.stats_p3p(k-1),g(:,:,k),g_estp(:,:,k),M*N);
            end
            obj.stats_p1.t_nmse=M*beta_dk(1)*sigma^2*K/(beta_dk(1)*P_C*T_1+sigma^2); %Not sure if T_1 goes on top, doesn't work without it, theory doesn't predict it
            obj.stats_p2.t_nmse=trace(pinv(P_C*Phi_2*pinv(psi_2)*Phi_2'+pinv(C_BI)));
            
            temp=0;
            for ii=1:(K-1)*ceil(N/M)
                phase_3(ii).k=ceil(ii/ceil(N/M))+1;
                phase_3(ii).first=(ii-floor(ii/ceil(N/M))*(ceil(N/M))-1)*M+1:(ii-floor(ii/ceil(N/M))*ceil(N/M)-1)*M+M;
                phase_3(ii).conditional=(floor(ii/ceil(N/M))~=ii/ceil(N/M));
                phase_3(ii).second=(ceil(N/M)-1)*M+1:N;
                if phase_3(ii).conditional==0
                    phase_3(ii).delta=phase_3(ii).second;
                else
                    phase_3(ii).delta= phase_3(ii).first;
                end
                for jj=1:size(phase_3(ii).delta,2)
                    phase_3(ii).Gp(:,jj)=g(:,phase_3(ii).delta(jj),1);
                    phase_3(ii).Lamb(jj)=obj.lamb(phase_3(ii).delta(jj),phase_3(ii).k);
                end
                phase_3(ii).Clamb=obj.C_lamb{ii};
                temp=trace(pinv(P_C*phase_3(ii).Gp'*pinv(psi_3)*phase_3(ii).Gp+pinv(phase_3(ii).Clamb)))+temp;
            end
            obj.stats_p3lp.t_nmse=temp/(K-1)*ceil(N/M);
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
            T_1=obj.T_1;
            K=obj.channel.Users;
            obj.s_1=hadamard(2^ceil(log2(T_1))); %generates orthogonal pilot signals
            for i=1:2^ceil(log2(T_1))-K %removes unecessary rows
                obj.s_1(2^ceil(log2(T_1))-i,:)=[];
            end
            obj.T_1=2^ceil(log2(T_1)); %number of training symbols in first phase
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
            obj.C_BI=C_BI;

            for ii =1:T_3
                %conditioning
                C_lamb{ii}=mean(temp2(ii).array,3);
                obj.C_lamb{ii}=C_lamb{ii};
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

            H_1=H_1+sqrt(beta_1/db2pow(rice_factor))*(randn(size(H_1)) + 1i*randn(size(H_1)))/sqrt(2);

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
            for k=2:K
                for ii=1:N
                    obj.lamb(ii,k,1)=g(1,ii,k)/g(1,ii,1);
                end
            end
        end

        function obj = setP1symbols(obj,T_1)
            obj.T_1=T_1;
        end

        function obj = setP2symbols(obj,T_2)
            obj.T_2=T_2;
        end

        function obj = setP3symbols(obj,T_3)
            obj.T_3=T_3;
        end
        
        function obj = channelCapacity(obj)
            H_0_est1 = obj.g_est(:,:,1);
            H_0_estk = obj.g_est(:,:,2);
            h_d_est1 = obj.h_dkest(:,1);
            h_d_estk = obj.h_dkest(:,2);
            channel = obj.channel;
            P = channel.power_Users;
            sigma = channel.noise_abs;

            v = exp(1j * phase(H_0_estk'*h_d_estk));
            obj.phi = v;
            obj = obj.find_effective_channel;
            H = obj.effective_channel(:,2);

            H_est = h_d_estk+H_0_estk*v;
            w_k = H_est/norm(H_est);

            SNR = P*abs(w_k'*(H))^2/sigma^2;
            user_rate_3P = log2(1+SNR)*(1-(obj.T_1+obj.T_2+obj.T_3)*channel.symbol_period/channel.coherence_period)/channel.symbol_period;

            v = zeros(size(H_0_estk,2),1);
            obj.phi = v;
            obj = obj.find_effective_channel;
            H = obj.effective_channel(:,1);

            H_dc = h_d_est1+H_0_est1*v;
            w_k = H_dc/norm(H_dc);

            SNR = P*abs(w_k'*(H))^2/sigma^2;
            user_rate_DC = log2(1+SNR)*(1-(obj.T_1)*channel.symbol_period/channel.coherence_period)/channel.symbol_period;
            
            obj.user_rate_3P = user_rate_3P;
            obj.user_rate_DC = user_rate_DC;
            %rate = (1-estimation_period/coherence_period)*capacity;
        end
    end
end

