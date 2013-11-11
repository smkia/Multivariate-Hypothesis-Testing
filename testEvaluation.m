function [sensitivity,specificity] = testEvaluation(result,mask)
TN = sum(sum(sum(and(~mask,~result)))); 
TP = sum(sum(sum(and(mask,result))));
sensitivity = TP / sum(sum(sum(mask)));
specificity = TN / sum(sum(sum(~mask)));
