%% [TIP7188 - Filtragem Adaptativa]
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