% ---------------------------------------------
% Roulette wheel selection algorithm
% input variables are selected candidate
% array and the probability array with
% same dimension. Only one candiate will
% be selected
% Example: a_i = Roulette_wheel(a,p)
% a = [2 4 5 8]; p = [0.2 0.7 0.05 0.05]
% a_i = 4 will be returned most frequently
% ---------------------------------------------
%
function a_selected = Roulette_wheel(a,p)
%
% dimension check, only 1D array are accpected
[row_a, col_a] = size(a);
[row_p, col_p] = size(p);
if (row_a~=1 && col_a~=1) || (row_p~=1 && col_p~=1)
    error('Only 1D arrays are accepted!');
elseif length(a) ~= length(p)
    error('Inconsistent dimensions of input arrays!');
end
%
accumulation = cumsum(p);
p_threshold = rand()*accumulation(end);
%
a_selected = a(accumulation > p_threshold);
if ~isempty(a_selected)
    a_selected = a_selected(1);
else
    a_selected = a(end);
end