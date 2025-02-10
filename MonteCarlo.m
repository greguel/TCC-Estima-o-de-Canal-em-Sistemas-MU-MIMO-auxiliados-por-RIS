classdef MonteCarlo
    properties
        n_reps
        n_par
        range
        mcStats_p1
        mcStats_p2
        mcStats_p3
        mcStats_p3p
        mcStats_p3l
        mcStats_p3lp
        mcStats_DFT1
        mcStats_DFT2
        mcStats_DFTLS1
        mcStats_DFTLS2
        channel
        mcStats
        mcNMSE_p1
        mcNMSE_p1t
        mcNMSE_p2
        mcNMSE_p2t
        mcNMSE_p3
        mcNMSE_p3l
        mcNMSE_p3lp
        mcNMSE_p3lpt
        mcNMSE_DFT1
        mcNMSE_DFT2
        mcNMSE_DFTLS1
        mcNMSE_DFTLS2
        mcNMSE_DFTt
        mcNMSE_pavg
        mcNMSE_p3p
        range_type
        method
        plot_avg_channel_3P
        plot_p3p_channel_3P
        plot_noise_SNR
        plot_LS_estimation
        plot_lambda_estimation
        DFTS
        DFTT_s
        P3T_1
        P3T_2
        P3T_3
        user_rate_DFT
        mcuser_rate_DFT
        user_rate_LS
        mcuser_rate_LS
        user_rate_DC
        mcuser_rate_DC
        mcuser_rate_3P
        user_rate_3P
        mcuser_rate_PE
        user_rate_PE
    end
    methods
        function obj=MonteCarlo(n_reps,channel)
            obj.n_reps=n_reps;
            obj.channel=channel;
            K=channel.Users;
            obj.mcStats_p1=repmat(Statistics, 1, n_reps);
            obj.mcStats_p2=repmat(Statistics, 1, n_reps);
            obj.mcStats_p3=cell(1,n_reps);
            obj.mcStats_p3p=cell(1,n_reps);
            obj.mcStats_p3l=repmat(Statistics, 1, n_reps);
            obj.mcStats_p3lp=repmat(Statistics, 1, n_reps);
            obj.mcStats_DFT1=repmat(Statistics, 1, n_reps);
            obj.mcStats_DFT2=repmat(Statistics, K, n_reps);
            obj.mcStats_DFTLS1=repmat(Statistics, 1, n_reps);
            obj.mcStats_DFTLS2=repmat(Statistics, K, n_reps);
            obj.mcStats=[];
            obj.mcNMSE_p1=[];
            obj.mcNMSE_p1t=[];
            obj.mcNMSE_p2=[];
            obj.mcNMSE_p2t=[];
            obj.mcNMSE_p3=[];
            obj.mcNMSE_DFT1=[];
            obj.mcNMSE_DFT2=[];
            obj.mcNMSE_DFTLS1=[];
            obj.mcNMSE_DFTLS2=[];
            obj.mcNMSE_DFTt=[];
            obj.mcNMSE_p3p=[];
            obj.mcNMSE_p3l=[];
            obj.mcNMSE_p3lp=[];
            obj.mcNMSE_p3lpt=[];
            obj.n_par=0;
            obj.range=[];
            obj.range_type=0;
            obj.mcNMSE_pavg=[];
            obj.method=0;
            obj.plot_avg_channel_3P=0;
            obj.plot_p3p_channel_3P=0;
            obj.plot_noise_SNR=0;
            obj.plot_LS_estimation=0;
            obj.plot_lambda_estimation=0;
            obj.DFTS=0;
            obj.DFTT_s=0;
            obj.P3T_1=0;
            obj.P3T_2=0;
            obj.P3T_3=0;
            obj.mcuser_rate_DFT=[];
            obj.user_rate_DFT=[];
            obj.mcuser_rate_LS=[];
            obj.user_rate_LS=[];
            obj.user_rate_DC=[]; %
            obj.mcuser_rate_DC=[];
            obj.mcuser_rate_3P=[];
            obj.user_rate_3P=[];
            obj.mcuser_rate_PE=[];
            obj.user_rate_PE=[];
        end
%rep functions
        function obj=mcRep3P(obj,i,simulation)
            simulation=simulation.makeChannels();
            simulation=simulation.firstPhase();
            simulation=simulation.secondPhase();
            simulation=simulation.thirdPhase();
            obj.user_rate_3P(i) = simulation.user_rate_3P;
            obj.user_rate_DC(i) = simulation.user_rate_DC; %
            simulation=simulation.simulationStats();
            obj.mcStats_p1(i)=simulation.stats_p1;
            obj.mcStats_p2(i)=simulation.stats_p2;
            obj.mcStats_p3l(i)=simulation.stats_p3l;
            obj.mcStats_p3lp(i)=simulation.stats_p3lp;
            obj.mcStats_p3{:,i}=simulation.stats_p3;
            obj.mcStats_p3p{:,i}=simulation.stats_p3p;
        end
        
        function obj=mcRepDFT(obj,i)
            simulation=SimulationDFT(obj.channel);
            simulation=simulation.makeChannels();
            if obj.DFTS~=0
                simulation=simulation.setSubphases(obj.DFTS);
            end
            if obj.DFTT_s~=0
                simulation=simulation.setNsymbols(obj.DFTT_s);
            end
            simulation=simulation.channelEstimation();
            obj.user_rate_DFT(i) = simulation.user_rate_DFT;
            obj.user_rate_LS(i) = simulation.user_rate_LS;
            obj.user_rate_PE(i) = simulation.user_rate_PE; %
            simulation=simulation.simulationStats();
            obj.mcStats_DFT1(i)=simulation.stats_DFT1;
            obj.mcStats_DFT2(:,i)=simulation.stats_DFT2;
            obj.mcStats_DFTLS1(i)=simulation.stats_DFTLS1;
            obj.mcStats_DFTLS2(:,i)=simulation.stats_DFTLS2;
        end

        function obj=mcNrep3P(obj)
            simulation=Simulation3P(obj.channel);
            if obj.P3T_1~=0
                simulation=simulation.setP1symbols(obj.P3T_1);
            end
            if obj.P3T_2~=0
                simulation=simulation.setP2symbols(obj.P3T_2);
            end
            if obj.P3T_3~=0
                simulation=simulation.setP3symbols(obj.P3T_3);
            end
            simulation=simulation.channel_avg;
            for i=1:obj.n_reps
                obj=obj.mcRep3P(i,simulation);
            end
        end
        
        function obj=mcNrepDFT(obj)
            for i=1:obj.n_reps
                obj=obj.mcRepDFT(i);
            end
        end

        %Change the channel
        function obj = mcProcess(obj, range, method, range_type)

            obj.range = range;
            n_par = size(obj.range, 2);
            obj.n_par = n_par;
            obj.mcStats = cell(5, n_par); % Assuming maximum 5 types of statistics

            obj.range_type=range_type;
            if strcmp(range_type, 'Distance')
                method_func = @changeUserdistance;
            elseif strcmp(range_type, 'RIS')
                method_func = @changeris;
            elseif strcmp(range_type, 'Users')
                method_func = @changeusers;
            elseif strcmp(range_type, 'BS')
                method_func = @changebs;
            elseif strcmp(range_type, 'Noise')
                method_func = @changenoise;
            end

            obj.method=method;
            for i = 1:n_par
                obj.channel = method_func(obj.channel, obj.range(i));
       
                switch method
                    case 'DFT'
                        obj = obj.mcNrepDFT();
                        obj.mcStats(4, i) = {obj.mcStats_DFT1};
                        obj.mcStats(5, i) = {obj.mcStats_DFT2};
                        obj.mcStats(7, i) = {obj.mcStats_DFTLS1};
                        obj.mcStats(8, i) = {obj.mcStats_DFTLS2};
                        obj.mcuser_rate_DFT(:,i)=obj.user_rate_DFT;
                        obj.mcuser_rate_LS(:,i)=obj.user_rate_LS;
                        obj.mcuser_rate_PE(:,i)=obj.user_rate_PE;
                        %

                    case '3P'
                        obj = obj.mcNrep3P();
                        obj.mcStats(1, i) = {obj.mcStats_p1};
                        obj.mcStats(2, i) = {obj.mcStats_p2};
                        obj.mcStats(3, i) = {obj.mcStats_p3};
                        obj.mcStats(6, i) = {obj.mcStats_p3p};
                        obj.mcStats(9, i) = {obj.mcStats_p3l};
                        obj.mcStats(10, i) = {obj.mcStats_p3lp};
                        obj.mcuser_rate_3P(:,i)=obj.user_rate_3P;
                        obj.mcuser_rate_DC(:,i)=obj.user_rate_DC;
                    case 'all'
                        obj = obj.mcNrep3P();
                        obj.mcStats(1, i) = {obj.mcStats_p1};
                        obj.mcStats(2, i) = {obj.mcStats_p2};
                        obj.mcStats(3, i) = {obj.mcStats_p3};
                        obj.mcStats(6, i) = {obj.mcStats_p3p};
                        obj.mcStats(9, i) = {obj.mcStats_p3l};
                        obj.mcStats(10, i) = {obj.mcStats_p3lp};
                        obj = obj.mcNrepDFT();
                        obj.mcStats(4, i) = {obj.mcStats_DFT1};
                        obj.mcStats(5, i) = {obj.mcStats_DFT2};
                        obj.mcStats(7, i) = {obj.mcStats_DFTLS1};
                        obj.mcStats(8, i) = {obj.mcStats_DFTLS2};
                        obj.mcuser_rate_DFT(:,i)=obj.user_rate_DFT;
                        obj.mcuser_rate_LS(:,i)=obj.user_rate_LS; %
                        obj.mcuser_rate_PE(:,i)=obj.user_rate_PE;
                        obj.mcuser_rate_3P(:,i)=obj.user_rate_3P;
                        obj.mcuser_rate_DC(:,i)=obj.user_rate_DC;
                end
                [i n_par]
            end

            switch method
                case 'DFT'
                    obj = obj.calcNMSEDFT();
                    obj = obj.calcUser_rate();
                case '3P'
                    obj = obj.calcNMSE3P();
                    obj = obj.calcUser_rate();
                case 'all'
                    obj = obj.calcNMSEDFT();
                    obj = obj.calcNMSE3P();
                    obj = obj.calcUser_rate();%
            end
        end
        %calcs NMSE
        function obj=calcUser_rate(obj)
            obj.mcuser_rate_DFT=mean(obj.mcuser_rate_DFT);
            obj.mcuser_rate_LS=mean(obj.mcuser_rate_LS);
            obj.mcuser_rate_PE=mean(obj.mcuser_rate_PE);
            obj.mcuser_rate_3P=mean(obj.mcuser_rate_3P);
            obj.mcuser_rate_DC=mean(obj.mcuser_rate_DC);
            %
        end

        function obj=calcNMSE3P(obj)
            mcStats=obj.mcStats;
            K=obj.channel.Users;
            n_reps=obj.n_reps;
            for i=1:size(mcStats,2)
                Stats_p1{i}=mcStats{1,i};
                Stats_p2{i}=mcStats{2,i};
                Stats_p3{i}=mcStats{3,i};
                Stats_p3p{i}=mcStats{6,i};
                Stats_p3l{i}=mcStats{9,i};
                Stats_p3lp{i}=mcStats{10,i};
            end
            for j=1:size(Stats_p1,2)
                for i=1:size(Stats_p1{j},2)
                    combined_nmsenum_p1{j}(:,:,i)=Stats_p1{j}(i).nmse_num;
                    combined_nmsedem_p1{j}(:,:,i)=Stats_p1{j}(i).nmse_dem;
                    combined_nmsenum_p1t{j}(i)=Stats_p1{j}(i).t_nmse;
                    combined_nmsenum_p2{j}(:,:,i)=Stats_p2{j}(i).nmse_num;
                    combined_nmsenum_p2t{j}(i)=Stats_p2{j}(i).t_nmse;
                    combined_nmsedem_p2{j}(:,:,i)=Stats_p2{j}(i).nmse_dem;
                    combined_nmsenum_p3l{j}(:,:,i)=Stats_p3l{j}(i).nmse_num;
                    combined_nmsenum_p3lp{j}(:,:,i)=Stats_p3lp{j}(i).nmse_num;
                    combined_nmsenum_p3lpt{j}(i)=Stats_p3lp{j}(i).t_nmse;
                    combined_nmsedem_p3l{j}(:,:,i)=Stats_p3l{j}(i).nmse_dem;
                    for k=1:size(Stats_p3{j}{i},2)
                        combined_nmsenum_p3{j}(:,:,k,i)=Stats_p3{j}{i}(k).nmse_num;
                        combined_nmsedem_p3{j}(:,:,k,i)=Stats_p3{j}{i}(k).nmse_dem;
                        combined_nmsenum_p3p{j}(:,:,k,i)=Stats_p3p{j}{i}(k).nmse_num;
                        combined_nmsedem_p3p{j}(:,:,k,i)=Stats_p3p{j}{i}(k).nmse_dem;
                    end
                end
            end
            for j=1:size(Stats_p1,2)
            obj.mcNMSE_p1(j)=trace(mean(combined_nmsenum_p1{j},3))/trace(mean(combined_nmsedem_p1{j},3));
            obj.mcNMSE_p1t(j)=mean(combined_nmsenum_p1t{j})/trace(mean(combined_nmsedem_p1{j},3));
            obj.mcNMSE_p2(j)=trace(mean(combined_nmsenum_p2{j},3))/trace(mean(combined_nmsedem_p2{j},3));
            obj.mcNMSE_p2t(j)=trace(mean(combined_nmsenum_p2t{j}))/trace(mean(combined_nmsedem_p2{j},3));
            obj.mcNMSE_p3l(j)=trace(mean(combined_nmsenum_p3l{j},3))/trace(mean(combined_nmsedem_p3l{j},3));
            obj.mcNMSE_p3lp(j)=trace(mean(combined_nmsenum_p3lp{j},3))/trace(mean(combined_nmsedem_p3l{j},3));
            obj.mcNMSE_p3lpt(j)=trace(mean(combined_nmsenum_p3lpt{j}))/trace(mean(combined_nmsedem_p3l{j},3));
            obj.mcNMSE_p3(j)=trace(mean(combined_nmsenum_p3{j},[3,4]))/trace(mean(combined_nmsedem_p3{j},[3,4]));
            obj.mcNMSE_p3p(j)=trace(mean(combined_nmsenum_p3p{j},[3,4]))/trace(mean(combined_nmsedem_p3p{j},[3,4]));
            re_combined_nmsenum_p2=reshape(combined_nmsenum_p2{j},[size(combined_nmsenum_p2{j},[1,2]),1,size(combined_nmsenum_p2{j},3)]);
            re_combined_nmsedem_p2=reshape(combined_nmsedem_p2{j},[size(combined_nmsedem_p2{j},[1,2]),1,size(combined_nmsedem_p2{j},3)]);
            obj.mcNMSE_pavg(j)=trace(mean(cat(3,re_combined_nmsenum_p2,combined_nmsenum_p3{j}),[3,4]))/trace(mean(cat(3,re_combined_nmsedem_p2,combined_nmsedem_p3{j}),[3,4]));
            end
        end


        function obj=calcNMSEDFT(obj)
            mcStats=obj.mcStats;
            for i=1:size(mcStats,2)
                Stats_DFT1{i}=mcStats{4,i};
                Stats_DFT2{i}=mcStats{5,i};
                Stats_DFTLS1{i}=mcStats{7,i};
                Stats_DFTLS2{i}=mcStats{8,i};
            end
            for j=1:size(Stats_DFT1,2)
                for i=1:size(Stats_DFT1{j},2)
                    combined_nmsenum_DFT1{j}(:,:,i)=Stats_DFT1{j}(i).nmse_num;
                    combined_nmsedem_DFT1{j}(:,:,i)=Stats_DFT1{j}(i).nmse_dem;
                    combined_nmsenum_DFTLS1{j}(:,:,i)=Stats_DFTLS1{j}(i).nmse_num;
                    combined_t_nmse_DFT1(i,j)=Stats_DFT1{j}(i).t_nmse;
                    for k=1:size(Stats_DFT2{j},1)
                        combined_nmsenum_DFT2{j}(:,:,k,i)=Stats_DFT2{j}(k,i).nmse_num;
                        combined_nmsenum_DFTLS2{j}(:,:,k,i)=Stats_DFTLS2{j}(k,i).nmse_num;
                        combined_nmsedem_DFT2{j}(:,:,k,i)=Stats_DFT2{j}(k,i).nmse_dem;
                        combined_t_nmse_DFT2(i,j,k)=Stats_DFT2{j}(1,i).t_nmse;
                    end
                end
            end
            for j=1:size(Stats_DFT1,2)
            obj.mcNMSE_DFT1(j)=trace(mean(combined_nmsenum_DFT1{j},3))/trace(mean(combined_nmsedem_DFT1{j},3));
            obj.mcNMSE_DFT2(j)=trace(mean(combined_nmsenum_DFT2{j},[3,4]))/trace(mean(combined_nmsedem_DFT2{j},[3,4]));
            obj.mcNMSE_DFTLS1(j)=trace(mean(combined_nmsenum_DFTLS1{j},3))/trace(mean(combined_nmsedem_DFT1{j},3));
            obj.mcNMSE_DFTLS2(j)=trace(mean(combined_nmsenum_DFTLS2{j},[3,4]))/trace(mean(combined_nmsedem_DFT2{j},[3,4]));
            obj.mcNMSE_DFTt(j,1)=mean(combined_t_nmse_DFT1(:,j),1);
            obj.mcNMSE_DFTt(j,2)=mean(combined_t_nmse_DFT2(:,j,:),[1,3]);
            end
        end

%plot
        function mcPlot(obj)
            channel=obj.channel;
            P_C=channel.power_Users;
            ttau=channel.symbol_period;
            beta_dk=channel.pathloss_BS_Users;
            range_type=obj.range_type;
            method=obj.method;
            plot_avg_channel_3P=obj.plot_avg_channel_3P;
            plot_p3p_channel_3P=obj.plot_p3p_channel_3P;
            plot_noise_SNR=obj.plot_noise_SNR;
            plot_LS_estimation=obj.plot_LS_estimation;
            plot_lambda_estimation=obj.plot_lambda_estimation;
            figure('Position', [100, 100, 800, 600]);
            %plot
            if (strcmp(range_type,'Noise') && plot_noise_SNR==1)
                x_axis_range=pow2db(P_C./db2pow(obj.range)); 
            elseif strcmp(range_type,'RIS')
                x_axis_range=obj.range.^2;
            else
                x_axis_range=obj.range;
            end

            switch method
                case 'DFT'
                    semilogy(x_axis_range,obj.mcNMSE_DFT1,'r*',x_axis_range,obj.mcNMSE_DFT2,'g*',x_axis_range,obj.mcNMSE_DFTt(:,1),'r',x_axis_range,obj.mcNMSE_DFTt(:,2),'g')
                     lgd =legend('Direct channel','BS-RIS-UE channel','Theoretical direct channel','Theoretical BS-RIS-UE channel');
                     plot_avg_channel_3P=0;
                     plot_p3p_channel_3P=0;
                     plot_lambda_estimation=0;
                     
                case '3P'
                    semilogy(x_axis_range,obj.mcNMSE_p1,'r*',x_axis_range,obj.mcNMSE_p2,'g*',x_axis_range,obj.mcNMSE_p1t,'r',x_axis_range,obj.mcNMSE_p2t,'g')
                    lgd =legend('Direct channel phase 1','BS-RIS-UE channel phase 2','Theoretical direct channel', 'Theoretical BS-RIS-UE channel phase 2');
                    currentLabels = legend().String;
                    hold on
                    semilogy(x_axis_range,obj.mcNMSE_p3,'b*')
                    currentLabels{end+1}='BS-RIS-UE channel phase 3';
                    lgd =legend(currentLabels);
                    plot_LS_estimation=0;
                case 'all'
                    semilogy(x_axis_range,obj.mcNMSE_DFT1,'r*',x_axis_range,obj.mcNMSE_DFT2,'g*',x_axis_range,obj.mcNMSE_DFTt(:,1),'r',x_axis_range,obj.mcNMSE_DFTt(:,2),'g')
                    hold on
                    semilogy(x_axis_range,obj.mcNMSE_p1,'r+',x_axis_range,obj.mcNMSE_p2,'g+',x_axis_range,obj.mcNMSE_p3,'b+')
                    lgd =legend('Direct channel MMSE-DFT','BS-RIS-UE channel MMSE-DFT','Theoretical direct channel MMSE-DFT','Theoretical BS-RIS-UE channel MMSE-DFT','Direct channel phase 1 3P','BS-RIS-UE channel phase 2 3P','BS-RIS-UE channel phase 3 3P');
            end
            %adds plot for flags
            hold on
            if  plot_avg_channel_3P==1
                currentLabels = legend().String;
                currentLabels{end+1}='Average BS-RIS-UE channel';
                Labelcolor='c';
                if strcmp(method,'3P')
                    Labelcolor=[Labelcolor,'*'];
                else
                    currentLabels{end}=[currentLabels{end}, ' 3P'];
                    Labelcolor=[Labelcolor,'+'];
                end
                semilogy(x_axis_range,obj.mcNMSE_pavg,Labelcolor)
                lgd =legend(currentLabels);
            end
            if  plot_p3p_channel_3P==1
                currentLabels = legend().String;
                currentLabels{end+1}='BS-RIS-UE channel phase 3 with perfect previous phases';
                Labelcolor='m';
                if strcmp(method,'3P')
                    Labelcolor=[Labelcolor,'*'];
                else
                    currentLabels{end}=[currentLabels{end}, ' 3P'];
                    Labelcolor=[Labelcolor,'+'];
                end
                semilogy(x_axis_range,obj.mcNMSE_p3p,Labelcolor)
                lgd =legend(currentLabels);
            end

            if plot_LS_estimation==1
                currentLabels = legend().String;
                currentLabels{end+1}='Direct channel LS-DFT';
                Labelcolor='mo';
                semilogy(x_axis_range,obj.mcNMSE_DFTLS1,Labelcolor , 'MarkerSize', 13)
                currentLabels{end+1}='BS-RIS-UE channel LS-DFT';
                Labelcolor='co';
                semilogy(x_axis_range,obj.mcNMSE_DFTLS2,Labelcolor)
                lgd =legend(currentLabels);

            end
            
            if  plot_lambda_estimation==1
                currentLabels = legend().String;
                currentLabels{end+1}='Lambda';
                Labelcolor='y';
                if strcmp(method,'3P')
                    Labelcolor=[Labelcolor,'*'];
                else
                    currentLabels{end}=[currentLabels{end}, ' 3P'];
                    Labelcolor=[Labelcolor,'+'];
                end
                semilogy(x_axis_range,obj.mcNMSE_p3l,Labelcolor)
                if plot_p3p_channel_3P==1
                    currentLabels{end+1}='Lambda with perfect previous phases';
                    Labelcolor='k';
                    if strcmp(method,'3P')
                        Labelcolor=[Labelcolor,'*'];
                    else
                        currentLabels{end}=[currentLabels{end}, ' 3P'];
                        Labelcolor=[Labelcolor,'+'];
                    end
                    semilogy(x_axis_range,obj.mcNMSE_p3lp,Labelcolor)
                    currentLabels{end+1}='Lambda Theoretical with perfect previous phases';
                    Labelcolor='k';
                    if strcmp(method,'3P')
                    else
                        currentLabels{end}=[currentLabels{end}, ' 3P'];
                    end
                    semilogy(x_axis_range,obj.mcNMSE_p3lpt,Labelcolor)
                end
                lgd = legend(currentLabels);
            end
            %finish plot
            ylabel('NMSE')
            title('RIS MIMO Channel Estimation NMSE')
            %for each range type
            switch range_type
                case 'Distance'
                    xlabel('User distance from origin (m)')
                    currentLabels = legend().String;


                    currentLabels{end+1}='RIS position';

                    RIS_pos=obj.channel.position_RIS;
                    yLimits = get(gca, 'YLim'); 
                    line([RIS_pos, RIS_pos], yLimits, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
                    
                    lgd =legend(currentLabels);
                case 'RIS'
                    xlabel('Number of RIS elements')
                case 'Users'
                    xlabel('Number of Users')
                case 'BS'
                    xlabel('Number of BS antenna')
                case 'Noise'
                    if plot_noise_SNR==0
                        xlabel('\sigma^2(dBW)')
                    else
                        xlabel('SNR(dB)')
                    end
            end
        end
        
        function mcPlot2(obj)
            channel=obj.channel;
            P_C=channel.power_Users;
            ttau=channel.symbol_period;
            beta_dk=channel.pathloss_BS_Users;
            range_type=obj.range_type;
            method=obj.method;
            plot_avg_channel_3P=obj.plot_avg_channel_3P;
            plot_p3p_channel_3P=obj.plot_p3p_channel_3P;
            plot_noise_SNR=obj.plot_noise_SNR;
            plot_LS_estimation=obj.plot_LS_estimation;
            plot_lambda_estimation=obj.plot_lambda_estimation;
            figure('Position', [100, 100, 900, 600]);
            %plot
            if (strcmp(range_type,'Noise') && plot_noise_SNR==1)
                x_axis_range=pow2db(P_C./db2pow(obj.range)); %not sure if symbol period should be taken into consideration
            elseif strcmp(range_type,'RIS')
                x_axis_range=obj.range.^2;
            else
                x_axis_range=obj.range;
            end
                    subplot(1,2,1)

                    mcNMSE_DFT1 = obj.mcNMSE_DFT1;
                    mcNMSE_DFT2 = obj.mcNMSE_DFT2;
                    mcNMSE_DFTt = obj.mcNMSE_DFTt;

                    semilogy(x_axis_range,mcNMSE_DFT1,'r*',x_axis_range,mcNMSE_DFT2,'go',x_axis_range,transpose(mcNMSE_DFTt(:,1)),'r',x_axis_range,transpose(mcNMSE_DFTt(:,2)),'g--')
                    yl = ylim;  
                    %ylim([1e-8, yl(2)]);
                    lgd1 = legend('Canal Direto MMSE-DFT','Canal concatenado MMSE-DFT','Canal direto teórico MMSE-DFT','Canal concatenado teórico MMSE-DFT');
                    title('Error de estimação normalizado para MMSE-DFT e LS-DFT');
                    set(lgd1, 'Location', 'southoutside', 'Orientation', 'vertical');
                    grid on
                    hold off
                    subplot(1,2,2)
                    mcNMSE_p1 = obj.mcNMSE_p1;
                    mcNMSE_p2 = obj.mcNMSE_p2;
                    mcNMSE_p3 = obj.mcNMSE_p3;

                    semilogy(x_axis_range,mcNMSE_p1,'r*',x_axis_range,mcNMSE_p2,'go',x_axis_range,mcNMSE_p3,'b+')
                    lgd2 =legend('Canal direto, fase 1 3P','Canal concatenado, fase 2 3P','Canal concatenado, fase 3 3P');
                    title('Error de estimação normalizado para 3P');
                    set(lgd2, 'Location', 'southoutside', 'Orientation', 'vertical');
                    xlabel('\sigma^2(dBW)')
                    yl = ylim; 
                    %ylim([1e-8, yl(2)]);
                    grid on
                    hold off
            %adds plot for flags
            if  (plot_avg_channel_3P==0)
                subplot(1,2,2)
                hold on
                currentLabels = legend().String;
                currentLabels{end+1}='Canal concatenado médio 3P';
                Labelcolor='c';
                Labelcolor=[Labelcolor,'s'];
                mcNMSE_pavg = obj.mcNMSE_pavg;
   
                semilogy(x_axis_range,mcNMSE_pavg,Labelcolor)
                lgd2 =legend(currentLabels);
                hold off
            end
            if  plot_p3p_channel_3P==0
                subplot(1,2,2)
                hold on
                currentLabels = legend().String;
                currentLabels{end+1}='Canal concatenado, fase 3, com fases anteriores "perfeitas" 3P';
                Labelcolor='m';
                Labelcolor=[Labelcolor,'d'];
                mcNMSE_p3p = obj.mcNMSE_p3p;

                semilogy(x_axis_range,mcNMSE_p3p,Labelcolor)
                lgd2 =legend(currentLabels);
                hold off
            end

            if plot_LS_estimation==1
                subplot(1,2,1)
                hold on
                currentLabels = legend().String;
                currentLabels{end+1}='Canal direto LS-DFT';
                Labelcolor='md';
                mcNMSE_DFTLS1 = obj.mcNMSE_DFTLS1;
                mcNMSE_DFTLS2 = obj.mcNMSE_DFTLS2;

                semilogy(x_axis_range,mcNMSE_DFTLS1,Labelcolor, 'MarkerSize',11)
                currentLabels{end+1}='Canal concatenado LS-DFT';
                Labelcolor='cs';
                semilogy(x_axis_range,mcNMSE_DFTLS2,Labelcolor)
                lgd1 =legend(currentLabels);
                hold off

            end
            
            if  plot_lambda_estimation==1
                currentLabels = legend().String;
                currentLabels{end+1}='Lambda';
                Labelcolor='y';
                if strcmp(method,'3P')
                    Labelcolor=[Labelcolor,'*'];
                else
                    currentLabels{end}=[currentLabels{end}, ' 3P'];
                    Labelcolor=[Labelcolor,'+'];
                end
                semilogy(x_axis_range,obj.mcNMSE_p3l,Labelcolor)
                if plot_p3p_channel_3P==1
                    currentLabels{end+1}='Lambda with perfect previous phases';
                    Labelcolor='k';
                    if strcmp(method,'3P')
                        Labelcolor=[Labelcolor,'*'];
                    else
                        currentLabels{end}=[currentLabels{end}, ' 3P'];
                        Labelcolor=[Labelcolor,'+'];
                    end
                    semilogy(x_axis_range,obj.mcNMSE_p3lp,Labelcolor)
                    currentLabels{end+1}='Lambda Theoretical with perfect previous phases';
                    Labelcolor='k';
                    if strcmp(method,'3P')
                    else
                        currentLabels{end}=[currentLabels{end}, ' 3P'];
                    end
                    semilogy(x_axis_range,obj.mcNMSE_p3lpt,Labelcolor)
                end
                lgd = legend(currentLabels);
            end
            %finish plot
            ylabel('NMSE')
            %title('RIS MIMO Channel Estimation NMSE')
            %for each range type
            switch range_type
                case 'Distance'
                    xlabel('Distância do usuário da origem (m)')
                    currentLabels = legend().String;

   
                    currentLabels{end+1}='Posição da RIS';

                    RIS_pos=obj.channel.position_RIS;
                    yLimits = get(gca, 'YLim');  % Get current y-axis limits
                    line([RIS_pos, RIS_pos], yLimits, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
                    
                    lgd =legend(currentLabels);
                    subplot(1,2,2)
                    xlabel('Distância do usuário da origem (m)')
                    currentLabels = legend().String;

                    currentLabels{end+1}='Posição da RIS';

                    RIS_pos=obj.channel.position_RIS;
                    yLimits = get(gca, 'YLim');  
                    line([RIS_pos, RIS_pos], yLimits, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
                    
                    lgd =legend(currentLabels);
                case 'RIS'
                    xlabel('Number of RIS elements')
                case 'Users'
                    xlabel('Number of Users')
                case 'BS'
                    xlabel('Number of BS antenna')
                case 'Noise'
                    if plot_noise_SNR==0
                        xlabel('\sigma^2(dBW)')
                    else
                        xlabel('SNR(dB)')
                    end
            end
        end
        function plot_pathloss(obj)
            range_type=obj.range_type;
            range=obj.range;
            channel=obj.channel;
            if strcmp(range_type, 'Distance')
                figure(2);
                for ii=1:size(range,2)
                    channel = changeUserdistance(channel, range(ii));
                    beta_d(ii)=channel.pathloss_BS_Users(1);
                    beta_c(ii)=channel.pathloss_BS_RIS*channel.pathloss_RIS_Users(1);
                end
                semilogy(range,beta_d,'r*--',range,beta_c,'go:')
                
                title("Perda de percurso em função da distância da ERB")
                hold on
               
                RIS_pos=obj.channel.position_RIS;
                yLimits = get(gca, 'YLim');  
                line([RIS_pos, RIS_pos], yLimits, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
                
                xlabel("Posição do usuário(m)")
                ylabel("Ganho do canal")
                lgd =legend("Perda de percurso do canal direto", "Perda de percurso do canal concatenado","Posição da RIS");
            else
                disp("Error: wrong range type")
            end
        end

        function obj = mcPlotCapacity(obj)
            range_type=obj.range_type;
            range=obj.range;
            DFT=obj.mcuser_rate_DFT;
            LS=obj.mcuser_rate_LS;
            PE=obj.mcuser_rate_PE;
            P3=obj.mcuser_rate_3P;
            DC=obj.mcuser_rate_DC;

            figure(3)
            plot(range,DFT, 'r*--', range, LS,'go-.', range, P3, 'bs:', range, DC, 'cd-')
            lgd =legend("MMSE-DFT", "LS-DFT", "3P usuário típico", "Estimação MMSE sem RIS");
            xlabel("User Position(m)")
            ylabel("Taxa de Dados(bit/s)")
            title("Taxa de dados de cada método")
            grid on

            switch range_type
                case 'Distance'
                    xlabel('Distância do usuário da origem (m)')
                    currentLabels = legend().String;
                    currentLabels{end+1}='Posição da RIS';

                    RIS_pos=obj.channel.position_RIS;
                    yLimits = get(gca, 'YLim');  
                    line([RIS_pos, RIS_pos], yLimits, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);

                    lgd =legend(currentLabels);
                case 'RIS'
                    xlabel('Number of RIS elements')
                case 'Users'
                    xlabel('Number of Users')
                case 'BS'
                    xlabel('Number of BS antenna')
                case 'Noise'
                        xlabel('\sigma^2(dBW)')
            end
        end

        function obj=flagPlotAvgCH3P(obj,flag)
            obj.plot_avg_channel_3P=flag;
        end


        function obj=flagPlotp3pCH3P(obj,flag)
            obj.plot_p3p_channel_3P=flag;
        end


        function obj=flagNoiseSNR(obj,flag)
            obj.plot_noise_SNR=flag;
        end

        function obj=flagLambdaEstimation(obj,flag)
            obj.plot_lambda_estimation=flag;
        end

        function obj=flagLSestimation(obj,flag)
            obj.plot_LS_estimation=flag;
        end

        function obj=fixSubphases(obj,S)
            obj.DFTS=S;
        end

        function obj=fixNsymbols(obj,T_s)
            obj.DFTT_s=T_s;
        end

        function obj=fixP1symbols(obj,T_1)
            obj.P3T_1=T_1;
        end

        function obj=fixP2symbols(obj,T_2)
            obj.P3T_2=T_2;
        end

        function obj=fixP3symbols(obj,T_3)
            obj.P3T_3=T_3;
        end

    end
end