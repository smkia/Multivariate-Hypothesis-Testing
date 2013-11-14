function [h] = FDR(pValues,alpha,method)
dim = size(pValues); 
pValues = pValues(:);
[pValuesSorted, sortedIndices] = sort(pValues);  
testNum = length(pValues);
switch method
    case 'BH'
        thresh = ((1:testNum)/testNum)  * alpha;
        thresh = thresh';
        h = (pValuesSorted<=thresh);  
    case 'BR'
        thresh = (((1:testNum)/testNum)/ c(testNum))  * alpha ;  
        thresh = thresh';
        h = (pValuesSorted<=thresh);  
    case 'BKY'
        alpha1 = alpha/(1+alpha);
        thresh = ((1:testNum)/testNum)  * alpha1;
        thresh = thresh';
        h = (pValuesSorted<=thresh);
        r1 = sum(h);
        if r1 ~= 0 || r1 ~= testNum
            alpha2 = (testNum/(testNum - r1)) * alpha1;
            thresh = ((1:testNum)/testNum)  * alpha2;
            thresh = thresh';
            h = (pValuesSorted<=thresh);
        end
end

% undo the sorting 
[~, unsort] = sort(sortedIndices); 
h = h(unsort);  
% convert the output back into the original format 
h = reshape(h, dim);
h = double(h);

function s = c(V) 
% See Genovese, Lazar and Holmes (2002) page 872, second column, first paragraph 
if V<1000   % compute it exactly   
    s = sum(1./(1:V)); 
else   % approximate it   
    s = log(V) + 0.57721566490153286060651209008240243104215933593992359880576723488486772677766467093694706329174674951463144724980708248096050401448654283622417399764492353625350033374293733773767394279259525824709491600873520394816567; 
end  