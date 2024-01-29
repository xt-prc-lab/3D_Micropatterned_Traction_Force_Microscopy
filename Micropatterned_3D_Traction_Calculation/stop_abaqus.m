%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                    Function used to stop abaqus' standard, in Linux, when the calculations are completed.                    %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       exp_type [string]: string that is apended at the beginning of the result file names.                                   %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%   Last Revison Date: 22/01/2024                                                                                              %
%   Author: Manuel Gomez Gonzalez                                                                                              %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = stop_abaqus(exp_type, varargin)

    % There is a bug in Abaqus for Linux where "standard" would not stop when the calculation is completed. The script 
    % stop_abaqus.m checks if the simulation is completed and kills "standard". With the following line, stop_abaqus.m is run 
    % in the background while the rest of the matlab script is run.

    if nargout
        varargout{nargout} = [];
    end

    wait_t = 1;         % Waiting time in seconds.

    % When "standard" finishes its calculations, it will create an *.sta file. Periodically look for this file and kill 
    % "standard" after it appears.
    while 1

        pause(wait_t*10);

        sta = dir([exp_type, '*.sta']);
        prt = dir([exp_type, '*.prt']);

        if (~isempty(sta)) && (~isempty(prt)) && (sta.datenum <= prt.datenum)

            pause(wait_t*10);
            system('killall -9 standard')
            delete([exp_type, '*.sta']);
            delete([exp_type, '*.prt']);
            datetime('now', 'TimeZone', 'local', 'Format', 'd-MMM-y HH:mm:ss')

        end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%