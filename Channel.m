classdef Channel
    properties
        BS_antennas %number of BS antenna
        RIS_elements %number of RIS elements
        RIS_columns %number of RIS columns
        RIS_rows %number of RIS rows
        Users %number of single antenna users
        speed_light=300000000
        frequency %operating frequency
        wave_length
        distance_BS_elements %distance between BS antennas (on grid)
        distance_RIS_elements %distance between RIS elements (on grid)
        position_RIS %x position of RIS
        height_BS %height of the BS
        height_RIS %RIS height
        distance_BS_RIS %assumes BS x position is 0
        position_Users %vector of user x coordenates, one value for each user
        distance_BS_Users %assumes user height is 0
        distance_RIS_Users %distance between RIS and each user
        elevation_BS_RIS %elevation angle between the RIS and BS
        power_Users
        coherence_period
        symbol_period
        noise_dB
        noise_abs
        RIS_correlation_matrix
        BS_correlation_matrix
        pathloss_exponent_LOS
        pathloss_exponent_NLOS
        reference_pathloss_RIS_Users;
        reference_pathloss_BS_Users
        reference_pathloss_BS_RIS
        pathloss_BS_RIS
        pathloss_BS_Users
        pathloss_RIS_Users
        BS_Users_propogation_loss
        rice_factor
    end
    methods
        %constructor
        function obj = Channel(varargin)
            if nargin==19 %all parameters
                obj.BS_antennas= varargin{1};
                obj.RIS_columns= varargin{2};
                obj.RIS_rows=varargin{3};
                obj.Users=varargin{4};
                obj.frequency=varargin{5};
                obj.position_RIS=varargin{6};
                obj.height_BS=varargin{7};
                obj.height_RIS=varargin{8};
                obj.position_Users=varargin{9};
                obj.power_Users=varargin{10};
                obj.coherence_period=varargin{11};
                obj.symbol_period=varargin{12};
                obj.noise_dB=varargin{13};
                obj.RIS_correlation_matrix=varargin{14};
                obj.BS_correlation_matrix=varargin{15};
                obj.pathloss_exponent_LOS=varargin{16};
                obj.pathloss_exponent_NLOS=varargin{17};
                obj.reference_pathloss_RIS_Users=varargin{18};
                obj.reference_pathloss_BS_RIS=varargin{19};
                obj.reference_pathloss_BS_Users=varargin{20};
                obj.BS_Users_propogation_loss=varargin{21};
                obj.rice_factor=varargin{22};

            end
            if nargin==6 %essential parameters with defaults
                obj.BS_antennas= varargin{1};
                obj.RIS_columns= varargin{2};
                obj.RIS_rows=varargin{3};
                obj.RIS_elements=obj.RIS_columns*obj.RIS_rows;
                obj.Users=varargin{4};
                obj.frequency=2.4*10^9;
                obj.position_RIS=50; 
                obj.height_BS=10;
                obj.height_RIS=20; 
                obj.position_Users=varargin{5};
                obj.power_Users=1;
                obj.coherence_period=0.1;
                obj.symbol_period=10*10^(-6);
                obj.noise_dB=varargin{6};
                obj.RIS_correlation_matrix=zeros(obj.RIS_elements,obj.RIS_elements,obj.Users);%RIS correlation matrix (NxNxK)
                obj.BS_correlation_matrix=zeros(obj.BS_antennas,obj.BS_antennas,obj.Users);%BS correlation matrix (MxMxK)
                for k=1:obj.Users
                    obj.RIS_correlation_matrix(:,:,k)=eye(obj.RIS_elements);%RIS correlation matrix (NxNxK)
                    obj.BS_correlation_matrix(:,:,k)=eye(obj.BS_antennas);%BS correlation matrix (MxMxK)
                end
                obj.pathloss_exponent_LOS=2.2;
                obj.pathloss_exponent_NLOS=3.67;
                obj.reference_pathloss_RIS_Users=25;
                obj.reference_pathloss_BS_Users=25;
                obj.reference_pathloss_BS_RIS=20;
                obj.BS_Users_propogation_loss=1; %attenuation due to NLOS
            end
            obj.RIS_elements=obj.RIS_columns*obj.RIS_rows;
            obj.wave_length=obj.speed_light/obj.frequency;
            obj.distance_BS_elements=obj.wave_length/2;
            obj.distance_RIS_elements=obj.wave_length/2;
            obj.distance_BS_RIS=sqrt(obj.position_RIS^2+(obj.height_BS-obj.height_RIS)^2);
            obj.distance_BS_Users=sqrt(obj.position_Users.^2+obj.height_BS^2);
            obj.distance_RIS_Users=sqrt((obj.position_Users-obj.position_RIS).^2+obj.height_RIS^2);
            obj.elevation_BS_RIS=atan(abs(obj.height_BS-obj.height_RIS)/obj.position_RIS);
            obj.noise_abs=sqrt(db2pow(obj.noise_dB));
            obj.pathloss_BS_RIS=10^((-obj.reference_pathloss_BS_RIS)/10)/(obj.distance_BS_RIS^obj.pathloss_exponent_LOS);
            obj.pathloss_BS_Users=10^((-obj.reference_pathloss_BS_Users)/10)./(obj.distance_BS_Users.^obj.pathloss_exponent_NLOS);
            obj.pathloss_RIS_Users=10^((-obj.reference_pathloss_RIS_Users)/10)./(obj.distance_RIS_Users.^obj.pathloss_exponent_NLOS);
            obj.rice_factor=7;
        end
       
        function obj = changenoise(obj,noise_dB)
            obj.noise_dB=noise_dB;
            obj.noise_abs=sqrt(db2pow(obj.noise_dB));
        end

        function obj = changebs(obj, n_BS)
            obj.BS_antennas=n_BS;
            obj.BS_correlation_matrix=zeros(obj.BS_antennas,obj.BS_antennas,obj.Users);
            for k=1:obj.Users
                obj.BS_correlation_matrix(:,:,k)=eye(obj.BS_antennas);%BS correlation matrix (MxMxK)
            end
        end

        function obj = changeusers(obj,n_users)
            obj.Users=n_users;
            obj.RIS_correlation_matrix=zeros(obj.RIS_elements,obj.RIS_elements,obj.Users);
            for k=1:obj.Users
                obj.RIS_correlation_matrix(:,:,k)=eye(obj.RIS_elements);%RIS correlation matrix (NxNxK)
                obj.BS_correlation_matrix(:,:,k)=eye(obj.BS_antennas);%BS correlation matrix (MxMxK)
            end
            obj.position_Users(n_users)=obj.position_Users(n_users-1);
            obj.distance_BS_Users=sqrt(obj.position_Users.^2+obj.height_BS^2);
            obj.distance_RIS_Users=sqrt((obj.position_Users-obj.position_RIS).^2+obj.height_RIS^2);
            obj.pathloss_BS_Users=10^((-obj.reference_pathloss_BS_Users)/10)./(obj.distance_BS_Users.^obj.pathloss_exponent_NLOS);
            obj.pathloss_RIS_Users=10^((-obj.reference_pathloss_RIS_Users)/10)./(obj.distance_RIS_Users.^obj.pathloss_exponent_NLOS);
        end

        function obj= changeris(obj,rissidesize)
            obj.RIS_columns= rissidesize;
            obj.RIS_rows=rissidesize;
            obj.RIS_elements=obj.RIS_columns*obj.RIS_rows;
            obj.RIS_correlation_matrix=zeros(obj.RIS_elements,obj.RIS_elements,obj.Users);
            for k=1:obj.Users
                obj.RIS_correlation_matrix(:,:,k)=eye(obj.RIS_elements);%RIS correlation matrix (NxNxK)
            end
        end

        function obj = changeUserdistance(obj, distance)
            obj.position_Users=distance*ones(1,obj.Users);
            obj.distance_BS_Users=sqrt(obj.position_Users.^2+obj.height_BS^2);
            obj.distance_RIS_Users=sqrt((obj.position_Users-obj.position_RIS).^2+obj.height_RIS^2);
            obj.pathloss_BS_Users=10^((-obj.reference_pathloss_BS_Users)/10)./(obj.distance_BS_Users.^obj.pathloss_exponent_NLOS)/obj.BS_Users_propogation_loss;
            obj.pathloss_RIS_Users=10^((-obj.reference_pathloss_RIS_Users)/10)./(obj.distance_RIS_Users.^obj.pathloss_exponent_NLOS);
        end

        function obj = setPropagation_loss(obj, propagation_loss_dB)
            obj.BS_Users_propogation_loss = db2pow(propagation_loss_dB);
            obj.pathloss_BS_Users=10^((-obj.reference_pathloss_BS_Users)/10)./(obj.distance_BS_Users.^obj.pathloss_exponent_NLOS)/obj.BS_Users_propogation_loss;
        end
    end
end