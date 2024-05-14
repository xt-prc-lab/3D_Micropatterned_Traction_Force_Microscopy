%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                         Calculate the location of the maximum of a 3D matrix, with subpixel accuracy.                        %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       im [3D matrix]: 3D matrix.                                                                                             %
%       sz [int scalar]: number of points around the maximum used to locate the peak with subpixel accuracy. I.e., a box of    %
%                        [2*sz+1, 2*sz+1, 2*sz+1] points around the maximum will be used.                                      %
%       method [int scalar]: For method=1, calculate the peak location, with subpixel accuracy, as the centroid of the         %
%                            [2*sz+1, 2*sz+1, 2*sz+1] box. The threshold is set automatically so that the external frame is    %
%                            always below the thhreshold.                                                                      %
%                            For method=2, use a gaussian fit to locate the peak with subpixel accuracy.                       %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       xc [doble scalar]: x coordinate (coordinate along the columns) of the peak, with sub-pixel accuracy.                   %
%       yc [doble scalar]: y coordinate (coordinate along the rows) of the peak, with sub-pixel accuracy.                      %
%       zc [doble scalar]: z coordinate (coordinate along the layers) of the peak, with sub-pixel accuracy.                    %
%                                                                                                                              %
%   Last Revison Date: 02/04/2024                                                                                              %
%   Based on x_cntr_2D.m, which calculates the centroid in 2D, by Xavier Trepat (2007).                                        %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%   References:                                                                                                                %
%       "Two-dimensional Gaussian regression for sub-pixel displacement estimation in particle image velocimetry or particle   %
%       position estimation in particle tracking velocimetry"; H. Nobach and M. Honkanen; Experiments in Fluids, March 2005,   %
%       Volume 38, pp 511–515. https://doi.org/10.1007/s00348-005-0942-3                                                       %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xc, yc, zc] = x_cntr_3D(im, sz, method)

    [Ny, Nx, Nz] = size(im);
    [max_val, max_ind] = max(im, [], 'all', 'omitnan', 'linear');

    [mj, mi, mk] = ind2sub([Ny, Nx, Nz], max_ind);

    if isnan(max_val)

        xc = NaN;
        yc = NaN;
        zc = NaN;

    else

        % Locate the peak as the 3D centroid.
        if method==1

            x0 = max(mi-sz-1,1);
            y0 = max(mj-sz-1,1);
            z0 = max(mk-sz-1,1);

            xf = min(mi+sz+1,Nx);
            yf = min(mj+sz+1,Ny);
            zf = min(mk+sz+1,Nz);

            im = im(y0:yf, x0:xf, z0:zf);

            if isempty(im)

                xc = NaN;
                yc = NaN;
                zc = NaN;

            elseif any(isnan(im), 'all') || any(isinf(im), 'all')

                % In case the 3D matrix is infinity or nan, use the integer value.
                xc = 0;
                yc = 0;
                zc = 0;

            else

                mask = true(size(im));
                mask(2:end-1, 2:end-1, 2:end-1) = 0;
                th1 = max(im(mask));
                im = im - th1;
                im(im<0) = 0;

                [grdi, grdj, grdk] = meshgrid(x0:xf, y0:yf, z0:zf);

                xc = sum(im(:).*grdi(:))/sum(im(:));
                yc = sum(im(:).*grdj(:))/sum(im(:));
                zc = sum(im(:).*grdk(:))/sum(im(:));

            end

            xc = -xc+Nx/2;
            yc = -yc+Ny/2;
            zc = -zc+Nz/2;

        end

        % Gaussian fit. Similar procedure as https://doi.org/10.1007/s00348-005-0942-3, but extended to use any number of 
        % points around the maximum.
        if method==2

            n0 = -sz;
            nf =  sz;
            m0 = -sz;
            mf =  sz;
            p0 = -sz;
            pf =  sz;

            % Check that the stencil does not surpass the available points. Recalculate n0, nf, m0, mf, p0 and pf accordingly.
            x0 = max(mi + n0, 1);
            y0 = max(mj + m0, 1);
            z0 = max(mk + p0, 1);

            xf = min(mi + nf, Nx);
            yf = min(mj + mf, Ny);
            zf = min(mk + pf, Nz);

            im = im(y0:yf, x0:xf, z0:zf);

            % Since we are normalizing the images, the correlation might be negative, and thus the logarithm becomes complex. 
            % In that case, use the centroid method. In case the correlation is infinity or nan, use the integer value.
            if isempty(im)

                xc = NaN;
                yc = NaN;
                zc = NaN;

            elseif any(isnan(im), 'all') || any(isinf(im), 'all')

                xc = 0;
                yc = 0;
                zc = 0;

            elseif any(im<=0, 'all')

                [xc, yc, zc] = x_cntr_3D(im, sz, 1);

            else

                n0 = x0 - mi;
                nf = xf - mi;
                m0 = y0 - mj;
                mf = yf - mj;
                p0 = z0 - mk;
                pf = zf - mk;

                im = log(im);

                [grdi, grdj, grdk] = meshgrid(n0:nf, m0:mf, p0:pf);

                s =   sum(im, 'all');
                si =  sum(grdi.*im, 'all');
                si2 = sum(grdi.*grdi.*im, 'all');
                sj =  sum(grdj.*im, 'all');
                sj2 = sum(grdj.*grdj.*im, 'all');
                sk =  sum(grdk.*im, 'all');
                sk2 = sum(grdk.*grdk.*im, 'all');
                sij = sum(grdi.*grdj.*im, 'all');
                sik = sum(grdi.*grdk.*im, 'all');
                sjk = sum(grdj.*grdk.*im, 'all');

                c001 = 6/((n0-nf-1)*(m0-mf-1)*(p0-pf)*(p0-pf-1)*(p0-pf-2))*( ...
                    (p0+pf)*( 3*(n0+nf)^2/( (n0-nf)*(n0-nf-2) ) + 3*(m0+mf)^2/( (m0-mf)*(m0-mf-2) ) + 5*(p0^2+4*p0*pf+pf^2+p0-pf)/( (p0-pf+1)*(p0-pf-3) ) + 45*(m0+mf)^2/( (m0-mf)*(m0-mf-2) )*(n0^2+4*n0*nf+nf^2+n0-nf)*(n0+nf)^2/( (n0-nf+1)*(n0-nf)*(n0-nf-2)*(n0-nf-3) ) + 1 )*s - ...
                    6*(n0+nf)/( (n0-nf)*(n0-nf-2) )*(p0+pf)*si - ...
                    6*(m0+mf)/( (m0-mf)*(m0-mf-2) )*(p0+pf)*sj - ...
                    2*( 3*(n0+nf)^2/( (n0-nf)*(n0-nf-2) ) + 3*(m0+mf)^2/( (m0-mf)*(m0-mf-2) ) + 15*(p0+pf)^2/( (p0-pf+1)*(p0-pf-3) ) + 1 )*sk + ...
                    30*(p0+pf)/( (p0-pf+1)*(p0-pf-3) )*sk2 - ...
                   12*(n0+nf)/( (n0-nf)*(n0-nf-2) )*(m0+mf)^3/( (m0-mf)^2*(m0-mf-2)^2 )*(p0+pf)*sij + ...
                    12*(n0+nf)/( (n0-nf)*(n0-nf-2) )*( 1 - (m0+mf)^2/( (m0-mf)*(m0-mf-2) )*(p0+pf)^2/( (p0-pf)*(p0-pf-2) ) )*sik + ...
                    12*(m0+mf)/( (m0-mf)*(m0-mf-2) )*sjk ...
                    ) ;

                c002 = -30/( (n0-nf-1)*(m0-mf-1)*(p0-pf+1)*(p0-pf)*(p0-pf-1)*(p0-pf-2)*(p0-pf-3) )*( ...
                    (p0^2+4*p0*pf+pf^2+p0-pf)*s - ...
                    6*(p0+pf)*sk + ...
                    6*sk2 ...
                    ) ;

                c010 = 6/( (n0-nf-1)*(m0-mf)*(m0-mf-1)*(m0-mf-2)*(p0-pf-1) )*( ...
                    (m0+mf)*( 3*(n0+nf)^2/( (n0-nf)*(n0-nf-2) ) + 5*(m0^2+4*m0*mf+mf^2+m0-mf)/( (m0-mf+1)*(m0-mf-3) ) + 3*(p0+pf)^2/( (p0-pf)*(p0-pf-2) )*( 15*(n0^2+4*n0*nf+nf^2+n0-nf)*(n0+nf)^2/( (n0-nf+1)*(n0-nf)*(n0-nf-2)*(n0-nf-3) )+1 ) +1 )*s - ...
                    6*(n0+nf)/( (n0-nf)*(n0-nf-2) )*(m0+mf)*si - ...
                    2*( 3*(n0+nf)^2/( (n0-nf)*(n0-nf-2) ) + 3*(p0+pf)^2/( (p0-pf)*(p0-pf-2) ) + 15*(m0+mf)^2/( (m0-mf+1)*(m0-mf-3) ) + 1 )*sj + ...
                    30*(m0+mf)/( (m0-mf+1)*(m0-mf-3) )*sj2 - ...
                    6*(m0+mf)*(p0+pf)/( (p0-pf)*(p0-pf-2) )*sk + ...
                    12*(n0+nf)/( (n0-nf)*(n0-nf-2) )*( 1 - (m0+mf)^2/( (m0-mf)*(m0-mf-2) )*(p0+pf)^2/( (p0-pf)*(p0-pf-2) ) )*sij + ...
                   12*(p0+pf)/( (p0-pf)*(p0-pf-2) )*sjk - ...
                    12*(n0+nf)/( (n0-nf)*(n0-nf-2) )*(m0+mf)*(p0+pf)^3/( (p0-pf)^2*(p0-pf-2)^2 )*sik ) ;

                c020 = -30/( (n0-nf-1)*(m0-mf+1)*(m0-mf)*(m0-mf-1)*(m0-mf-2)*(m0-mf-3)*(p0-pf-1) )*( ...
                    (m0^2+4*m0*mf+mf^2+m0-mf)*s - ...
                    6*(m0+mf)*sj + ...
                    6*sj2 ...
                    ) ;

                c100 = 2/( (n0-nf)*(n0-nf-1)*(n0-nf-2)*(m0-mf-1)*(p0-pf-1) ) * ( ...
                    3*(n0+nf)*( 15*( (n0+nf)^2-1 )/( (n0-nf+1)*(n0-nf-3) ) + 3*(m0+mf)^2/( (m0-mf)*(m0-mf-2) ) + 3*(p0+pf)^2/( (p0-pf)*(p0-pf-2) ) - 4 )*s - ...
                    6*( 15*(n0+nf)^2/( (n0-nf+1)*(n0-nf-3) ) + 3*(m0+mf)^2/( (m0-mf)*(m0-mf-2) ) + 3*(p0+pf)^2/( (p0-pf)*(p0-pf-2) )+1 )*si - ...
                    18*(n0+nf)*(m0+mf)/( (m0-mf)*(m0-mf-2) )*sj - ...
                    18*(n0+nf)*(p0+pf)/( (p0-pf)*(p0-pf-2) )*sk + ...
                    90*(n0+nf)/( (n0-nf+1)*(n0-nf-3) )*si2 + ...
                   32*(m0+mf)/( (m0-mf)*(m0-mf-2) )*sij + ...
                    32*(p0+pf)/( (p0-pf)*(p0-pf-2) )*sik ...
                    ) ;

                c200 = -30/( (n0-nf+1)*(n0-nf)*(n0-nf-1)*(n0-nf-2)*(m0-mf-1)*(n0-nf-3)*(p0-pf-1) )*( ...
                    (n0^2+4*n0*nf+nf^2+n0-nf)*s - ...
                    6*(n0+nf)*si + ...
                    6*si2 ...
                    ) ;

                c110 = -36/( (n0-nf)*(n0-nf-1)*(n0-nf-2)*(m0-mf)*(m0-mf-1)*(m0-mf-2)*(p0-pf-1) )*( ...
                    (n0+nf)*(m0+mf)*s - ...
                    2*(n0+nf)*sj - ...
                    2*(m0+mf)*si + ...
                    4*sij ...
                    ) ;

                c101 = -36/( (n0-nf)*(n0-nf-1)*(n0-nf-2)*(m0-mf-1)*(p0-pf)*(p0-pf-1)*(p0-pf-2) )*( ...
                   (n0+nf)*(p0+pf)*s - ...
                    2*(p0+pf)*si - ...
                    2*(n0+nf)*sk + ...
                    4*sik ...
                    ) ;

                c011 = -36/( (n0-nf-1)*(m0-mf)*(m0-mf-1)*(m0-mf-2)*(p0-pf)*(p0-pf-1)*(p0-pf-2) )*( ...
                    (m0+mf)*(p0+pf)*( 15*(n0^2+4*n0*nf+nf^2+n0-nf)*(n0+nf)^2/( (n0-nf+1)*(n0-nf)*(n0-nf-2)*(n0-nf-3) ) + 1 )*s - ...
                    2*(p0+pf)*sj - ...
                    2*(m0+mf)*sk - ...
                    4*(n0+nf)/( (n0-nf)*(n0-nf-2) )*(m0+mf)^2/( (m0-mf)*(m0-mf-2) )*(p0+pf)*sij - ...
                    4*(n0+nf)/( (n0-nf)*(n0-nf-2) )*(m0+mf)*(p0+pf)^2/( (p0-pf)*(p0-pf-2) )*sik + ...
                    4*sjk ...
                    ) ;

                xc = 1/2*(2*c001*c020*c101 - c001*c011*c110 - 4*c002*c020*c100 + 2*c002*c010*c110 + c011^2*c100 - c010*c011*c101)/ ...
                    (4*c002*c020*c200 - c002*c110^2 - c011^2*c200 + c011*c101*c110 - c020*c101^2);
                yc = 1/2*(2*c001*c011*c200 - c001*c101*c110 - 4*c002*c010*c200 + 2*c002*c100*c110 + c010*c101^2 - c011*c100*c101)/ ...
                    (4*c002*c020*c200 - c002*c110^2 - c011^2*c200 + c011*c101*c110 - c020*c101^2);
                zc = 1/2*(2*c010*c011*c200 - c010*c101*c110 - 4*c001*c020*c200 + 2*c020*c100*c101 + c001*c110^2 - c011*c100*c110)/ ...
                    (4*c002*c020*c200 - c002*c110^2 - c011^2*c200 + c011*c101*c110 - c020*c101^2);

            end

            xc = -(xc + mi) + Nx/2;
            yc = -(yc + mj) + Ny/2;
            zc = -(zc + mk) + Nz/2;

            if (abs(xc) >= Nx/2) || (abs(yc) >= Ny/2) || (abs(zc) >= Nz/2)

                [xc, yc, zc] = x_cntr_3D(im, sz, 1);

            end

        end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
