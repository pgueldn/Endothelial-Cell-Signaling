function newmap = custom_cmap(m)

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

bottom = [0 0.3 0.5];
botmiddle = [0 0.6 0.75];
middle = [1 1 1];
topmiddle = [1 0.4 0];
top = [0.75 0.25 0];

midpoint=0; 

lims = get(gca, 'CLim');
if (lims(1) < midpoint) && (lims(2) > midpoint)
    ratio = abs(lims(1)-midpoint) / (abs(lims(1)-midpoint) + lims(2)-midpoint);
    neglen = round(m*ratio);
    poslen = m - neglen;
    
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, neglen);
    newmap1 = zeros(neglen, 3);
    
    for i=1:3
        newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, poslen);
    newmap = zeros(poslen, 3);
    
    for i=1:3
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    newmap = [newmap1; newmap];
    
elseif lims(1) >= midpoint
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
else
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
end
