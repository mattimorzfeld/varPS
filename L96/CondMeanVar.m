
function [mxgy,sxgy] = CondMeanVar(xs,m,inds,C)

mx  = m(inds);
sxx = C(inds,inds);
% my  = muo([1:inds(1)-1  inds(end)+1:end]);
% syy  = C([1:inds(1)-1  inds(end)+1:end],...
%                 [1:inds(1)-1  inds(end)+1:end]);
% sxy = C(inds,[1:inds(1)-1 inds(end)+1:end]);
% mxgy = mx+sxy*(syy\(xs([1:inds(1)-1 inds(end)+1:end])-my));
% sxgy = sxx-sxy*(syy\sxy');



M= 5;
n = length(C);

if inds(1)-1-M<=0 && inds(end)+1+M <=n
    my  = m([1:inds(1)-1  inds(end)+1:inds(end)+1+M]);
    syy  = C([1:inds(1)-1  inds(end)+1:inds(end)+1+M],...
        [1:inds(1)-1  inds(end)+1:inds(end)+1+M]);
    sxy = C(inds,[1:inds(1)-1 inds(end)+1:inds(end)+1+M]);
    mxgy = mx+sxy*(syy\(xs([1:inds(1)-1 inds(end)+1:inds(end)+1+M])-my));
    sxgy = sxx-sxy*(syy\sxy');
elseif inds(1)-1-M<=0 && inds(end)+1+M >n
    my  = m([1:inds(1)-1  inds(end)+1:end]);
    syy  = C([1:inds(1)-1  inds(end)+1:end],...
        [1:inds(1)-1  inds(end)+1:end]);
    sxy = C(inds,[1:inds(1)-1 inds(end)+1:end]);
    mxgy = mx+sxy*(syy\(xs([1:inds(1)-1 inds(end)+1:end])-my));
    sxgy = sxx-sxy*(syy\sxy');
elseif inds(1)-1-M>0 && inds(end)+1+M>n
    my  = m([inds(1)-1-M:inds(1)-1 inds(end)+1:end]);
    syy  = C([inds(1)-1-M:inds(1)-1 inds(end)+1:end],...
        [inds(1)-1-M:inds(1)-1 inds(end)+1:end]);
    sxy = C(inds,[inds(1)-1-M:inds(1)-1 inds(end)+1:end]);
    mxgy = mx+sxy*(syy\(xs([inds(1)-1-M:inds(1)-1 inds(end)+1:end])-my));
    sxgy = sxx-sxy*(syy\sxy');
else
    my  = m([inds(1)-1-M:inds(1)-1 inds(end)+1:inds(end)+1+M]);
    syy  = C([inds(1)-1-M:inds(1)-1 inds(end)+1:inds(end)+1+M],...
        [inds(1)-1-M:inds(1)-1 inds(end)+1:inds(end)+1+M]);
    sxy = C(inds,[inds(1)-1-M:inds(1)-1 inds(end)+1:inds(end)+1+M]);
    mxgy = mx+sxy*(syy\(xs([inds(1)-1-M:inds(1)-1 inds(end)+1:inds(end)+1+M])-my));
    sxgy = sxx-sxy*(syy\sxy');
end