function L_s_midday = importfile(workbookFile,sheetName,startRow,endRow)
%IMPORTFILE Import data from a spreadsheet
%   L_s_midday = IMPORTFILE(FILE) reads data from the first worksheet in
%   the Microsoft Excel spreadsheet file named FILE and returns the data as
%   column vectors.
%
%   L_s_midday = IMPORTFILE(FILE,SHEET) reads from the specified worksheet.
%
%   L_s_midday = IMPORTFILE(FILE,SHEET,STARTROW,ENDROW) reads from the
%   specified worksheet for the specified row interval(s). Specify STARTROW
%   and ENDROW as a pair of scalars or vectors of matching size for
%   dis-contiguous row intervals. To read to the end of the file specify an
%   ENDROW of inf.
%
%	Non-numeric cells are replaced with: NaN
%
% Example:
%   L_s_midday = importfile('Sol_L_s_midday.xlsx','Sheet1',1,669);
%
%   See also XLSREAD.

% Auto-generated by MATLAB on 2018/03/15 16:40:47

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 3
    startRow = 1;
    endRow = 669;
end

%% Import the data
[~, ~, raw] = xlsread(workbookFile, sheetName, sprintf('B%d:B%d',startRow(1),endRow(1)));
for block=2:length(startRow)
    [~, ~, tmpRawBlock] = xlsread(workbookFile, sheetName, sprintf('B%d:B%d',startRow(block),endRow(block)));
    raw = [raw;tmpRawBlock]; %#ok<AGROW>
end
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
I = cellfun(@(x) ischar(x), raw);
raw(I) = {NaN};
data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
L_s_midday = data(:,1);

