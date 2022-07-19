%% MÉTODOS 
% [TIP7188 - Filtragem Adaptativa]
% Author: Lucas Abdalah
%
% filter_hw.m
% 
% filter_hw is a package developped for the Adaptative Filtering Course
% It is a way to make a compilation for all function
%
% CONTENT
% 
%   SAVE DATA TO TXT FILE 
%       filter_hw.MAT2TXT      - Write a matrix X into a txt file
%       filter_hw.TENSOR2TXT   - Write a 3D tensor X into a txt file
% 
%   PLACE HOLDER
% 

classdef filter_hw

methods(Static)


%% HOMEWORK 2 - PROBLEM 5
function hw2p5(varargin)
% FILTER_HW.HW2P5  Perfom the error surface propose on the Hw 2, problem 5
%
%
%   See also.

    if isempty(varargin)
        save_results = false;
    else
        save_results = varargin{1};
    end

    N = 25;
    w_lim = 100;
    w = [linspace(-w_lim,w_lim,N); linspace(-w_lim,w_lim,N)];
    [w_0, w_1] = meshgrid(w(1,:), w(2,:));
    J_surface = @(w_0, w_1) 24.40 - 4.*w_0 - 9.*w_1 + w_0.^2 + w_1.^2;
    J = J_surface(w_0, w_1);
    h = figure();
    surf(w_0, w_1, J, 'EdgeColor', 'none');
    colormap turbo;
    xlabel('$w_0$', 'FontSize', 16, 'interpreter', 'latex');
    ylabel('$w_1$', 'FontSize', 16, 'interpreter', 'latex');
    zlabel('$J$', 'FontSize', 16, 'interpreter', 'latex');
    view([-24.5036297640653 47.6514617014408]);
    colorbar('box', 'off');
    grid on;
    axis tight;
    filter_hw.export_fig(save_results, h, 'figures/hw2p5');
end

%% VERBOSE DETAILS
function export_fig(Activate, h, filename)
    if Activate
        savefig_tight(h, filename, 'both');
        filter_hw.verbose_save(filename);
    else
        pause(1)
        close(h);
    end
end

function verbose_save(filename)
    fprintf('Saving Results for:\n\t %s \n', filename);
end


%% SAVE DATA TO TXT FILE
function mat2txt(filename, X, permission, header)
% ND.MAT2TXT  Write a matrix X into a txt file
%   mat2txt(filename, X, 'w', header) - Overwite the file
%   mat2txt(filename, X, 'a', header) - Append to the file end
%
%   See also.
        [I, J] = size(X);
        fileID = fopen(filename, permission);
        fprintf(fileID, [repelem('-', strlength(header)+3), '\n', header, ...
                '\n', repelem('-', strlength(header)+3), '\n']);
        fprintf(fileID, 'X(%d, %d)\n', I, J);
            for ii = 1:I
                for jj = 1:J
                    fprintf(fileID, ' %2.0f', X(ii,jj));
                end
                fprintf(fileID, ';\n');
            end
        fprintf(fileID, '\n');
        fclose(fileID);
end



% end methods list
end

end