%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                       Calculate the 3D cross-correlation of two 3D matrices, through Fourier Transform.                      %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       a [3D matrix]: first 3D matrix to cross-correlate.                                                                     %
%       b [3D matrix]: second 3D matrix to cross-correlate.                                                                    %
%       pad ['str']: if 'yes', calculate the full size cross-correlation matrix, with size rounded to the next power of two,   %
%                    i. e. the size of the Fourier Transforms and correlation matrix is 2^nextpow2(size(a) + size(b)).         %
%                    Finally  cut the size of the cross-correlation matrix to size(a) + size(b).                               %
%                    Otherwise, don't expand them, and the cross-correlation matrix has the same size as a and b.              %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       c [3D matrix]: cross-correlation of a and b.                                                                           %
%                                                                                                                              %
%   Last Revison Date: 27/03/2024                                                                                              %
%   Based on xcorrf2.m, created by R. Johnson (Revision 1.0, Date 1995/11/27). The original xcorrf2.m has dissapeared from     %
%   the Mathworks repositories.                                                                                                %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A                                                                                                                    %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = xcorrf3(a,b,pad)

    if nargin == 1  % For autocorrelation.
        b = a;
        pad = 'no';
    elseif nargin==2
       pad='no';
    end

    [ma,na,ra] = size(a);
    [mb,nb,rb] = size(b);

    if size(a) ~= size(b)
        pad = 'yes';
    end

    % Make reverse conjugate of one array.
    b = conj(b(mb:-1:1,nb:-1:1,rb:-1:1));

    if strcmpi(pad,'yes')    % Use power of 2 transform lengths:

        mf = 2^nextpow2(ma+mb);
        nf = 2^nextpow2(na+nb);
        rf = 2^nextpow2(ra+rb);

        a = fftn(a, [mf, nf, rf]);
        b = fftn(b, [mf, nf, rf]);

        % In order to normalize the peak "close to" 1, we take into account that the dimensions of "a" and "b" have changed.
        denom = sqrt(sum(abs(a - mean(a, 'all')).^2, 'all')*sum(abs(b - mean(b, 'all')).^2, 'all'))/(mf*nf*rf);

        % Multiply transforms then inverse transform:
        c = real(ifftn(a.*b));

        %  Trim to standard size:
%         c = c(1:ma+mb-1, 1:na+nb-1, 1:rb+rb-1);
        c = c(1:ma+mb, 1:na+nb, 1:rb+rb);

        % Normalize:
        c = c/denom;

    else

        % denom = std(a(:))*std(b(:))*ma*na*ra ;
        % It is twice as fast to use my own implementation of the variance.
        denom = sqrt(sum((a - mean(a, 'all')).^2, 'all')*sum((b - mean(b, 'all')).^2, 'all'));

        c = ifftshift(real(ifftn(fftn(a).*fftn(b))))/denom;

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
