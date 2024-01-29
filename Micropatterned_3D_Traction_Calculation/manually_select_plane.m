%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%              Show a 3D view of the beads stack and allow the user to select the upper surface plane of the gel.              %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       File [structure]: File structure coming from the 3D PIV.                                                               %
%       Settings [structure]: Settings structure coming from the 3D PIV.                                                       %
%       im [matrix]: 3D stack to display.                                                                                      %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       n_plane [int]: index of the plane selected.                                                                            %
%                                                                                                                              %
%   Last Revison Date: 23/01/2024                                                                                              %
%   Author: Manuel Gomez Gonzalez                                                                                              %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function n_plane = manually_select_plane(File, Settings, im)

    % Grid of the image.
    xx = 1:((File.TractionSize.i-1)*(1 - Settings.Overlap) + 1)*Settings.Resolution;
    yy = 1:((File.TractionSize.j-1)*(1 - Settings.Overlap) + 1)*Settings.Resolution;
    zz = 1:((File.TractionSize.k-1)*(1 - Settings.Overlap_z) + 1)*Settings.Resolution_z;

    xlim_0 = 1;
    xlim_f = min(xlim_0 + Settings.Size_ROIx - 1, max(xx(:)));
    ylim_0 = 1;
    ylim_f = min(ylim_0 + Settings.Size_ROIy - 1, max(yy(:)));
    zlim_0 = 1;
    zlim_f = min(zlim_0 + Settings.num_planes - 1, max(zz(:)));

    % Image parameters:
    fig_scale = 2;

    % Normalize the image:
    im_min = prctile(im,  0.75, 'all');
    im_max = prctile(im, 99.25, 'all');
    im = (im - im_min)/(im_max - im_min);

    tit_font_size = 10*fig_scale;

    % Colormaps
    c_beads = colormap('gray');

    % Titles.
    tit = 'Choose the surface plane';

    % Create and scale the figure.
    fig1 = figure(1);
    clf(fig1);
    set(fig1, 'Visible', 'on');
    pos = get(fig1, 'Position');
    pos(3:4) = fig_scale*pos(4)*[1, 1];
    set(fig1, 'Position', pos);

    % Create three axes in the figure, plot the projections of the beads imamage on them.
    % Parameters of the axes distribution geometry in the figure.
    nz = size(im, 3);

    Dx = (xlim_f - xlim_0)*Settings.PixelSizeXY;
    Dy = (ylim_f - ylim_0)*Settings.PixelSizeXY;
    Dz = nz*Settings.PixelSizeZ;

    DL = 0.07;
    DC = DL/1.6;
    DS = DC/2;
    HS = DS/1.5;
    DSR = DL;
    DR = DS + HS + DSR;
    DB = DR;
    DSB = DSR;
%     DU = DL;

    T = (Dx + Dz)/(1-DL-DC-DR);

    % Create the three axes with the specific geometry.
    ax{1} = axes('Parent', fig1, 'Position', [DL        , DB+Dz/T+DC, Dx/T, Dy/T], 'Color', 'None', 'Box', 'off', ...
                 'XTick', '', 'Xcolor', 'None', 'YTick', '', 'Ycolor', 'None', 'YDir', 'Normal');
    ax{2} = axes('Parent', fig1, 'Position', [DL        , DB        , Dx/T, Dz/T], 'Color', 'None', 'Box', 'off', ...
                 'XTick', '', 'Xcolor', 'None', 'YTick', '', 'Ycolor', 'None', 'YDir', 'Reverse');
    ax{3} = axes('Parent', fig1, 'Position', [DL+Dx/T+DC, DB+Dz/T+DC, Dz/T, Dx/T], 'Color', 'None', 'Box', 'off', ...
                 'XTick', '', 'Xcolor', 'None', 'YTick', '', 'Ycolor', 'None', 'YDir', 'Normal');

    handles.fig = fig1 ;
    n_plane = round((zlim_f - zlim_0)/2);

    % Create the sliders allowing to modify the 3D cuts of the stack shown.
    handles.bx  = uicontrol('Parent', fig1, 'Style', 'slider', 'Units', 'Normalized', ...
                            'Position', [DL, DSB, Dx/T, HS]                              , 'value', (xlim_f - xlim_0)/2, ...
                            'min', xlim_0, 'max', xlim_f);

    handles.by  = uicontrol('Parent', fig1, 'Style', 'slider', 'Units', 'Normalized', ...
                            'Position', [DL+Dx/T+DC+Dz/T+DS, DSB+HS+DS+Dz/T+DC, HS, Dy/T], 'value', (ylim_f - ylim_0)/2, ...
                            'min', ylim_0, 'max', ylim_f);

    handles.bzx = uicontrol('Parent', fig1, 'Style', 'slider', 'Units', 'Normalized', ...
                            'Position', [DL+Dx/T+DS, DSB+HS+DS, HS, Dz/T]                , 'value', (zlim_f - zlim_0)/2, ...
                            'min', zlim_0, 'max', zlim_f);

    handles.bzy = uicontrol('Parent', fig1, 'Style', 'slider', 'Units', 'Normalized', ...
                            'Position', [DL+Dx/T+DC, DSB+HS+DS+Dz/T+(DC-DS-HS), Dz/T, HS], 'value', (zlim_f - zlim_0)/2, ...
                            'min', zlim_0, 'max', zlim_f);

    % Assign an action (call the function replot_slider) to the movement of the sliders.
    set(handles.bx , 'Callback', {@replot_slider, handles});
    set(handles.by , 'Callback', {@replot_slider, handles});
    set(handles.bzx, 'Callback', {@replot_slider, handles});
    set(handles.bzy, 'Callback', {@replot_slider, handles});

    handles.n_plane = round(handles.bzy.Value);

    % Set the default values to the UI data.
    guidata(handles.fig, handles);

%     bgcolor = fig1.Color ;

    sgtitle(fig1, tit, 'FontSize', tit_font_size);

    % Prevent the figure to resize and mess with its geometry.
    set(fig1, 'Resize', 'off');

    % Call the redrawing of all the panels of the figure with its default values.
    replot_slider([], [], handles);

    %--------------------------------------------------------------------------------------------------------------------------%
    %                  Re-draw all the panels and move all the sliders of the figure when one slider is moved.                 %
    %--------------------------------------------------------------------------------------------------------------------------%

    function replot_slider( slider, ev, handles )

        % handles.bzx and handles.bzy move together, but handles.bzx is inverted. When moving handles.bzx, the value of 
        % handles.bzy needs to be calculated.
        if slider == handles.bzx
            handles.bzy.Value = -handles.bzx.Value + handles.bzx.Min + handles.bzx.Max;
        end

        % Clear the xy-axes of the figure and re-plot the image plane and cross in it.
        cla(ax{1});
        h{1} = imagesc(xx*Settings.PixelSizeXY, yy*Settings.PixelSizeXY, ...
                       im(:, :, round(handles.bzy.Value)), 'Parent', ax{1});
        axis(ax{1}, 'image', 'tight');
        set(ax{1}, 'Color', 'None', 'Box', 'off', ...
            'XLim', [xlim_0, xlim_f]*Settings.PixelSizeXY, 'XTick', '', 'Xcolor', 'None', ...
            'YLim', [ylim_0, ylim_f]*Settings.PixelSizeXY, 'YTick', '', 'Ycolor', 'None', ...
            'YDir', 'Normal', 'CLim', [0, 1]);
        colormap(ax{1}, c_beads);
        hold(ax{1}, 'on');
        % Plot the cross indicating the position of the x and y sliders.
        plot(ax{1}, handles.bx.Value*[1, 1]*Settings.PixelSizeXY, [ylim_0, ylim_f]*Settings.PixelSizeXY, '--b');
        plot(ax{1}, [xlim_0, xlim_f]*Settings.PixelSizeXY, handles.by.Value*[1, 1]*Settings.PixelSizeXY, '--b');

        % Clear the xz-axes of the figure and re-plot the image plane and cross in it.
        cla(ax{2});
        h{2} = imagesc(xx*Settings.PixelSizeXY, zz*Settings.PixelSizeZ, ...
                       squeeze(im(:, round( handles.by.Value), :))', 'Parent', ax{2});
        axis(ax{2}, 'image', 'tight');
        set(ax{2}, 'Color', 'None', 'Box', 'off', ...
            'XLim', [xlim_0, xlim_f]*Settings.PixelSizeXY, 'XTick', '', 'Xcolor', 'None', ...
            'YLim', [zlim_0, zlim_f]*Settings.PixelSizeZ, 'YTick', '', 'Ycolor', 'None', ...
            'YDir', 'Reverse', 'CLim', [0, 1]);
        colormap(ax{2}, c_beads);
        hold(ax{2}, 'on');
        % Plot the cross indicating the position of the x and z sliders.
        plot(ax{2}, handles.bx.Value*[1, 1]*Settings.PixelSizeXY, [zlim_0, zlim_f]*Settings.PixelSizeZ, '--b');
        plot(ax{2}, [xlim_0, xlim_f]*Settings.PixelSizeXY, handles.bzy.Value*[1, 1]*Settings.PixelSizeZ, '--b');

        % Clear the yz-axes of the figure and re-plot the image plane and cross in it.
        cla(ax{3});
        h{3} = imagesc(zz*Settings.PixelSizeZ, yy*Settings.PixelSizeXY, ...
                       squeeze(im(round(handles.bx.Value), :, :)), 'Parent', ax{3});
        axis(ax{3}, 'image', 'tight');
        set(ax{3}, 'Color', 'None', 'Box', 'off', ...
            'XLim', [zlim_0, zlim_f]*Settings.PixelSizeZ, 'XTick', '', 'Xcolor', 'None', ...
            'YLim', [ylim_0, ylim_f]*Settings.PixelSizeXY, 'YTick', '', 'Ycolor', 'None', ...
            'YDir', 'Normal', 'CLim', [0, 1]);
        colormap(ax{3}, c_beads);
        hold(ax{3}, 'on');
        % Plot the cross indicating the position of the y and z sliders.
        plot(ax{3}, handles.bzy.Value*[1, 1]*Settings.PixelSizeZ, [ylim_0, ylim_f]*Settings.PixelSizeXY, '--b');
        plot(ax{3}, [zlim_0, zlim_f]*Settings.PixelSizeZ, handles.by.Value*[1, 1]*Settings.PixelSizeXY, '--b');

        drawnow;

        n_plane = round(handles.bzy.Value);

        guidata(handles.fig, handles);

    end

    %--------------------------------------------------------------------------------------------------------------------------%

    guidata( handles.fig, handles );

    % Block the execution of the code until the figure is closed.
    uiwait( fig1 );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%