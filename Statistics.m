classdef Statistics
    properties
        var
        avg
        nmse_num
        nmse_dem
        t_nmse
    end
    methods
        function obj=Statistics()
            obj.var=0;
            obj.avg=0;
            obj.nmse_num = 0;
            obj.nmse_dem = 0;
            obj.t_nmse=0;
        end


        function obj=stats(obj,actual,estimate,N) %runs stats
            obj=error_avg(obj,actual,estimate,N);
            obj=error_var(obj,actual,estimate,N);
            obj=calc_nmse_num(obj,actual,estimate);
            obj=calc_nmse_dem(obj,actual);
        end


        function obj=error_var(obj,actual,estimate,N) %calculates the variance
            K=size(size(actual),2);
            error=(abs(actual-estimate));
            for i=1:K
                error=sum(error);
            end
            obj.var=error^2/N;
        end


        function obj=error_avg(obj,actual,estimate,N) %calculates the average
            K=size(size(actual),2);
            error=(abs(actual-estimate));
            for i=1:K
                error=sum(error);
            end
            obj.avg=error/N;
        end


        function obj=calc_nmse_num(obj,actual,estimate)
            obj.nmse_num=(actual-estimate)*(actual-estimate)';
        end


        function obj=calc_nmse_dem(obj,actual)
             obj.nmse_dem=(actual)*(actual)';
        end

        
    end
end