function G = debundle(G, varargin)
%DEBUNDLE Divided edge bundling
%
% G = debundle(G)
%
% This function replicates the divided edge bundling algorithm, as
% described in: 
% 
% Selassie D, Heller B, Heer J (2011) Divided edge bundling for directional
% network data. IEEE Trans Vis Comput Graph 17:2354-2363 
%
% This version is pretty efficient by Matlab standards, but of course much
% slower than Selassie's own C code; be prepared to wait a bit for graphs
% with thousands of edges.
%
% Input variables:
%
%   G:              directed graph object, including G.Nodes.x and
%                   G.Nodes.y fields holding the coordinates of each node. 
%
% Optional input variables (passed as parameter/value pairs:
%
%   ks:             edge spring constant, controls stretchiness of edges
%                   (lower = stretchier) [0.5e-3] 
%
%   kC:             edge coulomb constant, controls how strongly edges are
%                   attracted to each other (higher = more attracted) [2e4] 
%
%   l:              edge lane width, for edges going in opposite directions
%                   [25] 
%
%   s:              edge coulomb decay, i.e. attractive force range [30.0]
%
%   f:              friction for velocity damping [0.8]
%
%   dt:             time step for first set of iterations [20]
%
%   compatability:  Type of compatibility metrics to apply (angle, scale,
%                   position, visibility, and connective compatability)
%                   'all':       all metrics
%                   'noconnect': all metrics except connective
%                   'none':      none
%
%   solver:         ode solver to use.  I played around with a few
%                   alternatives to the leapfrog Selassie used when I was
%                   first coding this, and I've left that code in, but
%                   really you should stick with the default; it's faster
%                   and more stable, and the only one fully tested.
%                   ['selassie'] 
%                   'selassie': leapfrog-ish solver with damping applied
%                               outside the ODE (dxdt = v; dvdt = F)
%                   'leapfrog': leapfrog solver with damping solved within
%                               the ODE (dxdt = v; dvdt = F - vf) 
%                   'ode45':    uses ode45 solver (dxdt = v; dvdt = F - vf)
%
%   edgefun:        function handle for transform function for edge
%                   weights, of the form y = fun(x), where x and y are both
%                   nedge x 1 vectors. If edge weights span several orders
%                   of magnitude, the default linear treatment of edge
%                   weights may not be ideal.  This can be used to apply a
%                   tranformation before the (0,1] normalization.  Default
%                   is @(x) x, i.e. no transformation.    
%
% Output variables:
%
%   G:              directed graph object, same as input but with three
%                   additional properties added to each Edge:
%                   x:            n x 1 array, where n is the final number
%                                 of control points used in bundling, x
%                                 coordinates of the edge path 
%                   y:            n x 1 array, y coordinates of the
%                                 edge path 
%                   BundleCompat: 1 x nedge vector, edge compatibility
%                                 value between this edge and each other
%                                 one (set to 0 for edges running in
%                                 opposite directions).

% Copyright 2015 Kelly Kearney

%---------------------------------------
% Parse input
%---------------------------------------

p = inputParser;
p.addParameter('ks', 0.5e-3, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('kC', 2.0e4,  @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('l',  25,     @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('s',  35.0,   @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('f',   0.8,   @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('dt', 20,     @(x) validateattributes(x, {'numeric'}, {'scalar'})); % Selassie uses 40, but my tests are unstable there

p.addParameter('solver',        'selassie',   @(x) validateattributes(x, {'char'}, {}));
p.addParameter('compatability', 'all',        @(x) validateattributes(x, {'char'}, {}));
p.addParameter('edgefun',       @(x) x,       @(x) validateattributes(x, {'function_handle'}, {}));
p.addParameter('plot',          false,        @(x) validateattributes(x, {'logical'}, {'scalar'}));

p.parse(varargin{:});
A = p.Results;

validateattributes(G, {'digraph'}, {}, 1);

A.solver = validatestring(A.solver, {'selassie', 'leapfrog', 'ode45'}, '', 'solver');
A.compatability = validatestring(A.compatability, {'all', 'noconnect', 'none'}, '', 'compatability');

% A.ks = 0.5e-3;
% A.kC = 2.0e4;
% A.l = 25;
% A.s = 35.0; % 30.0
% A.f = 0.8;  % 0.2
% A.solver = 'selassie';
% A.plot = false;
% A.dt = 20; % Selassie uses 40, but my tests are unstable there
% A.compatability = 'all';
% A.edgefun = [];
% 
% A = parsepv(A, varargin);

%---------------------------------------
% Setup
%---------------------------------------

[~,loc] = ismember(G.Edges.EndNodes, G.Nodes.Name);

% Remove any loops

isloop = loc(:,1) == loc(:,2);

if any(isloop)
    Gorig = G;
    Gloop = rmedge(G, find(~isloop));
    G = rmedge(G, find(isloop));
    [~,loc] = ismember(G.Edges.EndNodes, G.Nodes.Name);
end

nedge = numedges(G);
nnode = numnodes(G);

% Undirected unweighted version (for connectivity)

adj = adjacency(G);
U = graph(adj|adj'); 

% Scale charge for number of edges
 
A.kC = A.kC./sqrt(nedge); 

%---------------------------------------
% Preprocessing
%---------------------------------------

fprintf('Preprocessing...\n');
tic

xedge = G.Nodes.x(loc)';
yedge = G.Nodes.y(loc)';

% Scale edges to fit in 1000 x 1000 box

xlim = minmax(xedge);
ylim = minmax(yedge);
fac = 1000./(max(diff(xlim), diff(ylim)));

xedge = (xedge - xlim(1)).*fac;
yedge = (yedge - ylim(1)).*fac;

% Edges as vectors

pvec = [diff(xedge); diff(yedge)];
pnorm = sqrt(sum(pvec.^2));

xmid = mean(xedge,1);
ymid = mean(yedge,1);

% Compatibility measures

dmin  = distances(U);   % Node distance
dminedge = nan(nedge);  % Edge distance 

A.issame = zeros(nedge);
[A.ca, A.cc, A.cs, A.cp, A.cv] = deal(nan(nedge));

c = ConsoleProgressBar;
c.setMaximum(nedge);
c.start();

for ip = 1:nedge
    c.setValue(ip);
    c.setText(sprintf('%4d/%4d',ip, nedge));
    for iq = 1:ip

        % Same or opposite direction
        
        pdotq = pvec(:,ip)'*pvec(:,iq);
        
        A.issame(ip,iq) = pdotq > 0;
        
        % Connectivity compatibility
        
        ii = loc(ip,[1 2 1 2]);
        jj = loc(iq,[1 1 2 2]);
        idx = ii + (jj-1)*nnode;
        
        dminedge(ip,iq) = min(dmin(idx));
        
        % Angle compatibility
        
        A.ca(ip,iq) = abs(pdotq./(pnorm(ip)*pnorm(iq)));
        
        % Scale compatibility
        
        lavg = (pnorm(ip)+pnorm(iq))/2;
        A.cs(ip,iq) = 2./(lavg./min(pnorm(ip),pnorm(iq)) + max(pnorm(ip),pnorm(iq))./lavg);
        
        % Position compatiblity
        
        rmid = sqrt((xmid(ip)-xmid(iq))^2 + (ymid(ip)-ymid(iq)).^2);
        A.cp(ip,iq) = lavg./(lavg + rmid);
        
        % Visibility compatibility
        
        tohead = [xedge(1,ip)-xedge(1,iq); yedge(1,ip)-yedge(1,iq)];
        totail = [xedge(2,ip)-xedge(1,iq); yedge(2,ip)-yedge(1,iq)];
        nrm = pvec(:,iq)./pnorm(iq);
        
        headonother = nrm .* (nrm'*tohead);
        tailonother = nrm .* (nrm'*totail);
        I = [[xedge(1,iq); yedge(1,iq)] + headonother [xedge(1,iq); yedge(1,iq)] + tailonother];
        
        tohead = [xedge(1,iq)-xedge(1,ip); yedge(1,iq)-yedge(1,ip)];
        totail = [xedge(2,iq)-xedge(1,ip); yedge(2,iq)-yedge(1,ip)];
        nrm = pvec(:,ip)./pnorm(ip);
        
        headonother = nrm .* (nrm'*tohead);
        tailonother = nrm .* (nrm'*totail);
        
        J = [[xedge(1,ip); yedge(1,ip)] + headonother [xedge(1,ip); yedge(1,ip)] + tailonother];

        Imid = mean(I,2);
        Jmid = mean(J,2);
        Ilen = sqrt(sum(diff(I,[],2).^2));
        Jlen = sqrt(sum(diff(J,[],2).^2));
        if Ilen == 0 || Jlen == 0
            A.cv(ip,iq) = 0;
        else
            midQmidI = sqrt((xmid(iq)-Imid(1)).^2 + (ymid(iq)-Imid(2)).^2);
            midPmidJ = sqrt((xmid(ip)-Jmid(1)).^2 + (ymid(ip)-Jmid(2)).^2);
            A.cv(ip,iq) = min(max(0, 1 - 2*midQmidI/Ilen), max(0, 1 - 2*midPmidJ/Jlen));
        end

    end
end
A.cc = 1./(1 + dminedge);
c.stop();
fprintf('\n ');

switch A.compatability
    case 'all'
        A.compat = A.ca .* A.cc .* A.cs .* A.cp .* A.cv;
    case 'none'
        A.compat = ones(nedge);
    case 'noconnect'
        A.compat = A.ca .* A.cs .* A.cp .* A.cv;
end

% Reflect lower triangle to fill symmetric matrices

A.issame = A.issame | A.issame';
A.compat = nansum(cat(3, A.compat, A.compat'), 3) ./ (eye(nedge)+1);

% Normalize edge weight to [0 1]

if isempty(A.edgefun)
    A.weight = G.Edges.Weight./max(G.Edges.Weight);
else
    A.weight = A.edgefun(G.Edges.Weight);
    A.weight = A.weight./max(A.weight);
end

toc

%---------------------------------------
% Bundling
%---------------------------------------

fprintf('Bundling...\n');


% Setup for integration

xc = xedge;
yc = yedge;

dt = A.dt;
nstep = 30;

% Integrate forward, using the 5-pass method ("Magic Iteration")
% decribed in Section 4.2 of Salassie et al. 2011.  

for ii = 1:5
    tic
    fprintf(' Pass %d\n',ii);

    % Split each edge segment in two

    xmid = (xc(1:end-1,:) + xc(2:end,:))./2;
    ymid = (yc(1:end-1,:) + yc(2:end,:))./2;

    nc = size(xc,1)+size(xmid,1);
    [xcnew, ycnew] = deal(zeros(nc,nedge));
    xcnew(1:2:end,:) = xc;
    xcnew(2:2:end,:) = xmid;
    ycnew(1:2:end,:) = yc;
    ycnew(2:2:end,:) = ymid;
    xc = xcnew;
    yc = ycnew;

    switch A.solver
        case 'selassie' % (simpler leapfrog... damping outside of equations)
            sz = size(xc);

            [x,y,vx,vy,ax,ay] = deal(nan([sz nstep]));
            x(:,:,1)  = xc;
            y(:,:,1)  = yc;
            vx(:,:,1) = 0; 
            vy(:,:,1) = 0;
            ax(:,:,1) = 0;
            ay(:,:,1) = 0;

            c = ConsoleProgressBar;
            c.setMaximum(nstep-1);
            c.start();

            for jj = 1:nstep-1
                c.setValue(jj);
                c.setText(sprintf('Pass %d: %2d of %2d', ii, jj, nstep-1));

                vx(:,:,jj+1) = (vx(:,:,jj) + ax(:,:,jj) .* dt/2) .* A.f;
                vy(:,:,jj+1) = (vy(:,:,jj) + ay(:,:,jj) .* dt/2) .* A.f;

                x(:,:,jj+1) = x(:,:,jj) + vx(:,:,jj+1).*dt;
                y(:,:,jj+1) = y(:,:,jj) + vy(:,:,jj+1).*dt;

                F = debforces(x(:,:,jj+1), y(:,:,jj+1), A, nedge, nc);
                ax(:,:,jj+1) = F(:,:,1);
                ay(:,:,jj+1) = F(:,:,2);

                vx(:,:,jj+1) = vx(:,:,jj+1) + ax(:,:,jj+1) .* dt/2;
                vy(:,:,jj+1) = vy(:,:,jj+1) + ay(:,:,jj+1) .* dt/2;
            end
            c.stop();
            fprintf('\n');

            xc = x(:,:,end);
            yc = y(:,:,end);

            if A.plot
                hfig = figure;
                h = plot(x(:,:,1), y(:,:,1), '-o');
                axis equal;
                for it = 1:nstep
                    for ie = 1:nedge
                        h(ie).XData = x(:,ie,it);
                        h(ie).YData = y(:,ie,it);
                    end
                    title(sprintf('Pass %d: Step %d of %d', ii, it, nstep));
                    pause
                end
                close(hfig);                    

            end

        case 'leapfrog' % dxdt = v; dvdt = F - vf 

            sz = size(xc);
            ts = (0:nstep-1)*dt;
            [x,y,vx,vy] = deal(nan([sz nstep]));
            x(:,:,1)  = xc;
            y(:,:,1)  = yc;
            vx(:,:,1) = 0; % index jj = jj+1/2
            vy(:,:,1) = 0;

            for jj = 1:nstep-1

                fac1 = 1 + A.f*dt/2;
                fac2 = (1 - A.f*dt/2).*(1 + A.f*dt/2);

                F = debforces(x(:,:,jj), y(:,:,jj), A, nedge, nc);

                vx(:,:,jj+1) = F(:,:,1).*dt.*fac1 + vx(:,:,jj).*fac2;
                vy(:,:,jj+1) = F(:,:,2).*dt.*fac1 + vy(:,:,jj).*fac2;

                x(:,:,jj+1) = x(:,:,jj) + vx(:,:,jj+1).*dt;
                y(:,:,jj+1) = y(:,:,jj) + vy(:,:,jj+1).*dt;

            end  

            xc = x(:,:,end);
            yc = y(:,:,end);

            % Debugging animation

            if A.plot
                hfig = figure;
                h = plot(x(:,:,1), y(:,:,1), '-o');
                axis equal;
                for it = 1:nstep
                    for ie = 1:nedge
                        h(ie).XData = x(:,ie,it);
                        h(ie).YData = y(:,ie,it);
                    end
                    title(sprintf('Pass %d: Step %d of %d', ii, it, nstep));
                    pause
                end
                close(hfig);
            end

        case 'ode45'

            xy0 = cat(3, xc, yc, zeros(size(xc)), zeros(size(yc)));
            xy0 = xy0(:);

            ts = [0 (nstep-1)*dt];

            [tout, yout] = ode45(@(t,xy) derivatives(t,xy,nc,nedge,A), ts, xy0);
            nt = size(yout,1);
            dtavg = mean(diff(tout));

            yout = reshape(yout, nt, nc, nedge, 4);
            xc = permute(yout(end,:,:,1), [2 3 1]);
            yc = permute(yout(end,:,:,2), [2 3 1]);
    end

    % Half the time step

    dt = dt/2;

    fprintf(' ');
    toc
end


%---------------------------------------
% Postprocessing
%---------------------------------------

fprintf('Postprocessing...\n');

% G.Edges.x = num2cell(xc, 1)';
% G.Edges.y = num2cell(yc, 1)';
% G.Edges.BundleCompat = num2cell(A.issame & A.compat > 0.05, 2);  

% Unscale

xc = xc./fac + xlim(1);
yc = yc./fac + ylim(1);

% for ii = 1:nedge
%     G.Edges.x{ii} = G.Edges.x{ii}/fac + xlim(1);
%     G.Edges.y{ii} = G.Edges.y{ii}/fac + ylim(1);
% end

% Add back any self-edges as repeating points

if any(isloop)
    [~, lloc] = ismember(Gloop.Edges.EndNodes(:,1), G.Nodes.Name);
    xcl = repmat(Gloop.Nodes.x(lloc)', nc, 1);
    ycl = repmat(Gloop.Nodes.y(lloc)', nc, 1);
    
    xc = [xc xcl];
    yc = [yc ycl];
    
    n = size(xc,2);
    tmp = A.compat;
    tmp(~A.issame) = 0;
    bcompat = zeros(n);
    bcompat(1:nedge,1:nedge) = tmp;
    
    srctar = [G.Edges.EndNodes; Gloop.Edges.EndNodes];
    order = findedge(Gorig, srctar(:,1), srctar(:,2));
    
    xc(:,order) = xc;
    yc(:,order) = yc;
    bcompat(order,order) = bcompat;
    
    G = Gorig;

else
    bcompat = A.compat;
    bcompat(~A.issame) = 0;
end
   
% Assign to graph object

G.Edges.x = num2cell(xc, 1)';
G.Edges.y = num2cell(yc, 1)';
G.Edges.BundleCompat = num2cell(bcompat, 2);  

fprintf('Done\n');


%---------------------------------------
% ODE function for 2D, 2nd order system
%---------------------------------------

% The variables in this system include x and y position and velocity for
% each control point on each edge. 

function da = derivatives(t, a, nc, nedge, A)

% Unpack the inputs

a = reshape(a, nc, nedge, 4); % control points x edge x coordinate

xc  = a(:,:,1);
yc  = a(:,:,2);
vxc = a(:,:,3);
vyc = a(:,:,4);

% Calculate forces

f = debforces(xc, yc, A, nedge, nc);

% ODE values

dxdt = vxc;                 % set dx/dt = x component of velocity
dydt = vyc;                 % set dy/dt = y component of velocity

dvxdt = f(:,:,1) - vxc*A.f; % dx2/d2t = x component of summed forces minus damping
dvydt = f(:,:,2) - vyc*A.f; % dy2/d2t = y component of summed forces minus damping

da = [dxdt(:); dydt(:); dvxdt(:); dvydt(:)];  % return the derivatives

%---------------------------------------
% Control point forces
%---------------------------------------

function F = debforces(xc, yc, A, nedge, nc)

F = zeros(nc,nedge,2);

% Calculate offset potential minimum locations 

xy = permute(cat(3, xc,yc), [3 2 1]); % xy x edge x control points
T = xy(:,:,3:end) - xy(:,:,1:end-2); % vector from previous control pint to following one
Tnorm = sqrt(sum(T.^2,1));

mj = zeros(2,nc,nedge);
for ic = 2:nc-1
    for iq = 1:nedge
        N = [0 -1; 1 0]*(T(:,iq,ic-1)./Tnorm(:,iq,ic-1)); %   norm(T(:,iq,ic-1)));
        mj(:,ic,iq) = [xc(ic,iq); yc(ic,iq)] + A.l*N;
    end
end

% Spring forces from neighbors

v1 = xy(:,:,1:end-2) - xy(:,:,2:end-1); % Vector from preceeding to current
v2 = xy(:,:,3:end  ) - xy(:,:,2:end-1); % Vector from next to current

r1 = sqrt(sum(v1.^2,1));
r2 = sqrt(sum(v2.^2,1));

v1unit = bsxfun(@rdivide, v1, r1);
v2unit = bsxfun(@rdivide, v2, r2);

Fs1 = bsxfun(@times, bsxfun(@times, A.ks .* nc .* r1, v1unit), A.weight');
Fs2 = bsxfun(@times, bsxfun(@times, A.ks .* nc .* r2, v2unit), A.weight');
Fs = Fs1 + Fs2;

% Add coulomb force from other corresponding control points

for ie = 1:nedge
    for ic = 2:nc-1
        
        p  = [xc(ic,  ie); yc(ic,  ie)]; % p_i
         
        mjtmp = xy(:,:,ic);                                          % Same direction
        mjtmp(:, ~A.issame(ie,:)) = mj(:,nc-(ic-1),~A.issame(ie,:)); % Opposite
        
        v = bsxfun(@minus, p, mjtmp);
        r = sqrt(sum(v.^2,1));
        vunit = bsxfun(@rdivide, v, r);
        
        Fcr =  -A.s.*A.kC .* r ./ (pi .* nc .* (A.s^2 + r.^2).^2);
        Fc = bsxfun(@times, vunit, Fcr);
        Fc(:,r == 0) = 0;
        
        Fc = bsxfun(@times, Fc, A.compat(ie,:).*A.weight');
        
        % Total
        
        F(ic,ie,:) = Fs(:,ie,ic-1) + sum(Fc,2);
        
    end
end



