function [h, Data] = plotdeb(G, varargin)
%PLOTDEB Plot divided edge bundled graph edges
%
% h = plotdeb(G, p1, v1, ...)
% [h, Data] = plotdeb(G, p1, v1, ...)
%
% Plot a graph using color-transitioning patches for edges; primarily
% intended for visualizing graphs following application of the divided edge
% bundling algorithm.
%
% Input variables:
%
%   G:          digraph object, including x and y path coordinates for
%               edges (or x and y coordinates for nodes when ititial mode
%               is used)  
%
% Optional input arguments (passed as parameter/value pairs):
%   
%   w:          width to apply to edges in initial mode.  Value assumes
%               data will be scaled to a 1000 x 1000 box, to maintain
%               consistency with Selassie's code. [10]  
%
%   p:          edge width fall-off [1.00]
%
%   wmin:       minimum edge width [0]
%
%   ax:         axis to plot to [gca]
%
%   smooth:     logical scalar, true to apply Gaussian smoothing to the
%               edge line paths [false]
%
%   alpha:      alpha value of gaussian filter used to smooth edge paths.
%               As alpha increases, the width of the window decreases (0 =
%               moving average) [1] 
%
%   initial:    logical scalar, true to plot straight-line, unbundled edges
%               [false]
%
%   selfloop:   Show edges that share the same node for source and target
%               as loops [true]
%
%   rloop:      Radius of circles used for self-loop edges [20]
%
%   gmax:       edge weight corresponding to maximum width value.  If
%               NaN, the maximum edge width in the graph will be used.
%               [NaN]
%
%   cthresh:    compatability threshold value over which two edges will
%               show visibly-summed weight [0.05]
%
% Output variables:
%
%   h:          graphics handle to the patch object
%
%   Data:       structre holding x, y, and c coordinates for the patch
%               object.

% Copyright 2015 Kelly Kearney

%--------------------
% Setup
%--------------------

validateattributes(G, {'digraph'}, {});

nedge = numedges(G);

p = inputParser;
p.addParameter('w',         10, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('p',          1, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('wmin',       0, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('alpha',      1, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('rloop',     20, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('gmax',     NaN, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('cthresh', 0.05, @(x) validateattributes(x, {'numeric'}, {'scalar'}));

p.addParameter('smooth',  false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('initial', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
p.addParameter('selfloop', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));

p.addParameter('ax', []);
p.addParameter('upsample', []);

p.parse(varargin{:});
Opt = p.Results;

% 
% Opt.w = 10;
% Opt.p = 1;
% Opt.wmin = 0;
% Opt.ax = [];
% Opt.smooth = false;
% Opt.alpha = 1;
% Opt.initial = false;
% Opt.selfloop = true;
% Opt.rloop = 20;
% Opt.gmax = NaN;
% % Opt.plot = true;
% % Opt.edgefun = [];
% Opt.upsample = [];
% 
% Opt = parsepv(Opt, varargin);

if ~Opt.initial && ~all(ismember({'x','y'}, G.Edges.Properties.VariableNames))
    Opt.initial = true;
end

if Opt.initial
    Opt.smooth = false;
end

%--------------------
% Calculate edge 
% details
%--------------------

% Scale factor for this dataset

xlim = minmax(G.Nodes.x);
ylim = minmax(G.Nodes.y);
fac = 1000./(max(diff(xlim), diff(ylim)));

% If initial, use node coordinates for edge paths

if Opt.initial
    nidx = reshape(findnode(G, G.Edges.EndNodes), size(G.Edges.EndNodes));
    
    x = G.Nodes.x(nidx)';
    y = G.Nodes.y(nidx)';
    
else
    x = cell2mat(G.Edges.x');
    y = cell2mat(G.Edges.y');
end

% Smooth edges if required

if Opt.smooth
    
    gfilt = gausswin(7, Opt.alpha);
    gfilt = gfilt./sum(gfilt);
    
    xtmp = filter(gfilt, 1, x, nan(6,1));
    ytmp = filter(gfilt, 1, y, nan(6,1));
    
    x(4:end-3,:) = xtmp(7:end,:);
    y(4:end-3,:) = ytmp(7:end,:);    
    
end

if ~isempty(Opt.upsample)
    [xnew, ynew] = deal(nan(Opt.upsample, size(x,2)));
    t = linspace(0,1,Opt.upsample);
    for ie = 1:size(x,2)
        xy = unique([x(:,ie), y(:,ie)], 'rows', 'stable');
        if size(xy,1) > 1
            pt = interparc(t, xy(:,1), xy(:,2), 'linear');
            xnew(:,ie) = pt(:,1);
            ynew(:,ie) = pt(:,2);
        else
            xnew(:,ie) = xy(1,1);
            ynew(:,ie) = xy(1,2);
        end
    end
    x = xnew;
    y = ynew;
end

% Bundle weight and visible edge thickness

if isnan(Opt.gmax)
    Opt.gmax = max(G.Edges.Weight);
end

gedge = G.Edges.Weight./Opt.gmax; % Normalized edge weight 

if Opt.initial
    gcontrol = ones(2,1) * gedge'; % Bundle weight of control points
else    
    gcontrol = zeros(size(x));
    
    d1 = (Opt.w.*gedge.^Opt.p)./fac; % Visual weight without bundling
    
    mask = cell2mat(G.Edges.BundleCompat) >= Opt.cthresh & ~eye(nedge);
    npt = size(x,1);
    for ie = 1:nedge
        
        xcmp = x(:,mask(ie,:));
        ycmp = y(:,mask(ie,:));
        gcmp = ones(npt,1) * gedge(mask(ie,:))';
        
        %---
        
        debug = false;
        
        if debug
            
            h1 = plot(xcmp, ycmp, '.');
            hold on;
            h2 = plot(x(:,ie), y(:,ie), 'marker', 'x', 'markersize', 10);
            arrayfun(@(x,y) rectangle('position', [x-d1(ie)/2 y-d1(ie)/2 d1(ie) d1(ie)], 'curvature', [1 1]), x(:,ie), y(:,ie));
          
        end
        
        D = ipdm([x(:,ie) y(:,ie)], [xcmp(:) ycmp(:)], ...
            'Result', 'Structure', 'Subset', 'Maximum', 'Limit', d1(ie)/2);

        gcontrol(:,ie) = gedge(ie);
        if ~isempty(D.rowindex)
            [~,cidx] = ind2sub(size(xcmp), D.columnindex);
            
            [ridx,tmp] = aggregate(D.rowindex, [cidx D.columnindex],  @(x) unique(x, 'rows'));
            wadd = nan(size(ridx));
            for ii = 1:length(tmp)
                [~,iunq] = unique(tmp{ii}(:,1));
                wadd(ii) = sum(gcmp(tmp{ii}(iunq,2)));
            end
           
            gcontrol(ridx,ie) = gcontrol(ridx,ie) + wadd;
        end
        
        if debug
            h3 = plot(xcmp(D.columnindex), ycmp(D.columnindex), 'o');
            dtmp = (Opt.w.*gcontrol(:,ie).^Opt.p)./fac;
            h4 = arrayfun(@(x,y,d) rectangle('position', [x-d/2 y-d/2 d d], 'curvature', [1 1]), x(:,ie), y(:,ie), dtmp, 'uni', 0);
            h4 = cat(1, h4{:});
            set(h4, 'edgecolor', 'r');
            
            gtmp = ones(npt,1)*gedge(ie);
            [tmpridx,gaddtmp] = aggregate(D.rowindex, gcmp(D.columnindex), @unique);
            gaddtmp = catuneven(2, 0, gaddtmp{:})';
            gadd = zeros(npt,size(gaddtmp,2));
            gadd(tmpridx,:) = gaddtmp;
            
%             figure;
%             bar([gtmp gadd], 'stacked');
            
        end
        
%         gcontrol = gcontrol./Opt.gmax; % Bundle weight of control points
    end
end


% Add self-loops (which aren't included in bundle weights, so just reflect
% their own weight regardless of overlap)

if Opt.selfloop
    nth = 21;
    th = linspace(pi, 3*pi, nth);
    
    isloop = strcmp(G.Edges.EndNodes(:,1), G.Edges.EndNodes(:,2));

    xloop = Opt.rloop/fac .* cos(th)' + Opt.rloop/fac;
    yloop = Opt.rloop/fac .* sin(th)';
    
    [xl,yl,gl] = deal(nan(nth,nedge));
    
    xl(:,isloop) = bsxfun(@plus, x(1,isloop), xloop);
    yl(:,isloop) = bsxfun(@plus, y(1,isloop), yloop);
    gl(:,isloop) = ones(size(xloop))*gedge(isloop)';
end

% Visual weight

d  = max((Opt.w .* gcontrol.^Opt.p)./fac, Opt.wmin./fac);
dl = max((Opt.w .* gl.^Opt.p)./fac, Opt.wmin./fac);

% Convert to polygons

% if Opt.selfloop
%     th = linspace(pi, 3*pi, 13);
%     xloop = Opt.rloop/fac .* cos(th)' + Opt.rloop/fac;
%     yloop = Opt.rloop/fac .* sin(th)';
% end
% isloop = strcmp(G.Edges.EndNodes(:,1), G.Edges.EndNodes(:,2));

[xp,yp,cp,wp] = deal(cell(nedge,1));
for ii = 1:nedge
    if Opt.selfloop && isloop(ii)
        [xp{ii},yp{ii},cp{ii},wp{ii}] = line2poly(xl(:,ii), yl(:,ii), dl(:,ii), gl(:,ii)*Opt.gmax); 
    else
        [xp{ii},yp{ii},cp{ii},wp{ii}] = line2poly(x(:,ii), y(:,ii), d(:,ii), gcontrol(:,ii)*Opt.gmax); 
    end
end

%--------------------
% Plot
%--------------------

if isempty(Opt.ax)
    Opt.ax = gca;
end

hold(Opt.ax, 'on');

[xp,yp,cp,wp] = singlepatch(xp,yp,cp,wp);

h = patch(xp,yp,cp,'parent', Opt.ax);
set(h, 'FaceAlpha', 0.4, 'edgecolor', 'none');   

Data.x = xp;
Data.y = yp;
Data.c = cp;
Data.w = wp;

% h = cellfun(@(x,y,c) patch(x,y,c, 'parent', Opt.ax), x,y,c, 'uni', 0);
% h = cat(1, h{:});

%--------------------
% Subfunction: 
% line2poly
%--------------------

function [xp,yp,cp,wp] = line2poly(x, y, w, g)
% x = x coordinate
% y = y coordinate
% w = width
% g = weight being depicted by that width

nx = length(x);
p = [x y]';
w = w/2;

% Unit normal for each line segment

v = diff(p, [], 2);
N = zeros(2,nx-1);
for iv = 1:nx-1
    N(:,iv) = [0 -1; 1 0]*(v(:,iv)./norm(v(:,iv))); % Unit normal right at each vector
end

% Check direction at each turn

c = zeros(3,nx);
for iv = 2:nx-1
    c(:,iv) = cross([p(:,iv+1) - p(:,iv-1); 0], [p(:,iv)-p(:,iv-1); 0]);
end

% Distance along line at each point

dr = sqrt(v(1,:).^2 + v(2,:).^2);
dr = [0 cumsum(dr)];
dr = dr./max(dr);

% Calculate new points (o = out segment, i = in segment)

lo = p(:,1:end-1) + bsxfun(@times, N, w(1:end-1)');  
ro = p(:,1:end-1) - bsxfun(@times, N, w(1:end-1)');  
li = p(:,2:end) + bsxfun(@times, N, w(2:end)');    
ri = p(:,2:end) - bsxfun(@times, N, w(2:end)');

li = [nan(2,1) li];
ri = [nan(2,1) ri];
lo = [lo nan(2,1)];
ro = [ro nan(2,1)];

% Calculate where intersections/mitering needs to be added

[xl,yl,xr,yr,cl, cr] = deal(cell(nx,1));
xl{1} = lo(1,1);
yl{1} = lo(2,1);
xr{1} = ro(1,1);
yr{1} = ro(2,1);
cl{1} = dr(1);
cr{1} = dr(1);

xl{end} = li(1,end);
yl{end} = li(2,end);
xr{end} = ri(1,end);
yr{end} = ri(2,end);
cl{end} = dr(end);
cr{end} = dr(end);

for ii = 2:nx-1
    
    if c(3,ii) < 1e-10    % Straight, for all intents and purposes
        xr{ii} = ri(1,ii);
        yr{ii} = ri(2,ii);
        xl{ii} = li(1,ii);
        yl{ii} = li(2,ii);
        cr{ii} = dr(ii);
        cl{ii} = dr(ii);
    else
        
        x1 = [lo(1,ii-1) li(1,ii) ri(1,ii) ro(1,ii-1) lo(1,ii-1)];
        y1 = [lo(2,ii-1) li(2,ii) ri(2,ii) ro(2,ii-1) lo(2,ii-1)];
        x2 = [lo(1,ii) li(1,ii+1) ri(1,ii+1) ro(1,ii) lo(1,ii)];
        y2 = [lo(2,ii) li(2,ii+1) ri(2,ii+1) ro(2,ii) lo(2,ii)];

        [xint, yint] = polyxpoly(x1, y1, x2, y2);
        if length(xint) == 1
            warning('Issue with polygon intersection... check this');
        end
        [~,imax] = max((xint - x(ii)).^2 + (yint - y(ii)).^2);
    
        if c(3,ii) < 0 % Turns left
       
            xl{ii} = xint(imax);
            yl{ii} = yint(imax);

            xr{ii} = [ri(1,ii); ro(1,ii)];
            yr{ii} = [ri(2,ii); ro(2,ii)];

            cl{ii} = dr(ii);
            cr{ii} = [dr(ii); dr(ii)];

        else           % Turns right

            xr{ii} = xint(imax);
            yr{ii} = yint(imax);

            xl{ii} = [li(1,ii); lo(1,ii)];
            yl{ii} = [li(2,ii); lo(2,ii)];

            cr{ii} = dr(ii);
            cl{ii} = [dr(ii); dr(ii)];

        end
    end
end

% Concatenate left and right

xl = cat(1, xl{:});
yl = cat(1, yl{:});
xr = cat(1, xr{:});
yr = cat(1, yr{:});
cl = cat(1, cl{:});
cr = cat(1, cr{:});

xp = [xl; xr(end:-1:1); xl(1)];
yp = [yl; yr(end:-1:1); yl(1)];
cp = [cl; cr(end:-1:1); cl(1)];

wp = interp1(dr', g, cp);

function visedgewidth()

g = G.Edges.Weight./max(G.Edges.Weight); 

if Opt.initial
    gnew = ones(2,1) * g';
else    
    gnew = zeros(size(x));
    
    d1 = (Opt.w.*g.^Opt.p)./fac; % Visual weight without bundling
    
    mask = cell2mat(G.Edges.BundleCompat) & ~eye(nedge);
    npt = size(x,1);
    for ie = 1:nedge
        
        xcmp = x(:,mask(ie,:));
        ycmp = y(:,mask(ie,:));
        gcmp = ones(npt,1) * g(mask(ie,:))';

        D = ipdm([x(:,ie) y(:,ie)], [xcmp(:) ycmp(:)], ...
            'Result', 'Structure', 'Subset', 'Maximum', 'Limit', d1(ie));

        gnew(:,ie) = g(ie);
        if ~isempty(D.rowindex)
            [~,cidx] = ind2sub(size(xcmp), D.columnindex);
            [ridx,b] = aggregate(D.rowindex, cidx, @unique);
            wadd = cellfun(@(x) sum(gcmp(x)), b);
            gnew(ridx,ie) = gnew(ridx,ie) + wadd;
        end
       
    end
end

gnew = gnew./max(gnew(:)); 
d = max((Opt.w .* gnew.^Opt.p)./fac, Opt.wmin./fac);







