%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                         Function used to modify the number of parallel workers of the local cluster.                         %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       n_w [optional, scalar]: number of workers.                                                                             %
%       n_th_per_worker [optional, scalar]: number of threads that each worker can use.                                        %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       [optional, scalar]: 1 if the function was successful, 0 otherwise.                                                     %
%                                                                                                                              %
%   Last Revison Date: 22/01/2024                                                                                              %
%   Author: Manuel Gomez Gonzalez                                                                                              %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = modify_workers_local_cluster(varargin)   % n_w, n_th_per_worker)

    try

        % Check the number of inputs.
        if nargin==0
            n_w = feature('NumCores');
            n_th_per_worker = n_w;

        elseif nargin==1
            n_w = varargin{1};
            n_th_per_worker = n_w;

        else
            n_w = varargin{1};
            n_th_per_worker = varargin{2};

        end

        % Change the maximum numbe of workers of the local cluster and the parellel pool.
        myCluster = parcluster('local') ;
        poolobj = gcp('nocreate') ;

        if myCluster.NumWorkers ~= n_w
            myCluster.NumWorkers = n_w;
            saveProfile(myCluster);
        end

        if isempty(poolobj) || poolobj.NumWorkers ~= n_w

            % Delete any pre-existing pool session.
            delete(poolobj);

            % Create the pool object.
            parpool(n_w);

        end

        pctRunOnAll(['maxNumCompThreads(', num2str(n_th_per_worker), ') ;']);

        % Return 1 if everything was successful.
        if nargout>=1
            varargout{1} = 1;
            varargout{nargout} = [];
        end

    catch

        % Return 0 if there was any problem.
        if nargout>=1
            varargout{1} = 0;
            varargout{nargout} = [];
        end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%