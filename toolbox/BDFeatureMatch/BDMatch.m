%--------------------------------------------------------------------------
% Bounded Distortion Feature Matching
%
% BDMatch.m
% This function filter candidate matches via the Bounded Distortion
% Alogirthm.
% Input: Porg (2xN) - each column has a pair of x,y coords from first image
%        Qorg (2xN) - coordinates of corresponding points in second image
%        K - ditostion bound of the deformation
%
% optionals: useboundary, displaydebug, createTemplateOptimizer
%
%
%
% This code, written by Eli Libman, advised by Meirav Galun,
% is an implementation of the algorithm described in the paper
% Feature Matching with Bounded Distortion
% by Yaron Lipman, Stav Yagev, Roi Poranne, David W. Jacobs, Ronen Basri
% ACM Trans. Graph., Volume 33 Issue 3, May 2014
%
% Disclaimer: The code is provided as-is and without any guarantees.
%--------------------------------------------------------------------------


function [ Pf, Qf, W, deformation, debugFigure ] = BDMatch( Porg,Qorg, K , pnorm, varargin )
%BDMATCH Bounded-Distortion correspondece filtering
%   This function filter candidate matches via the Bounded Distortion
% Alogirthm Input: P - 2*N - each column has a pair of x,y coords from
% first image Q - 2*N - coordinates of corresponding points in second image
% distortionBound - K bound constraining the deformation
%
% optionals: useboundary, displaydebug, createTemplateOptimizer
yalmip('clear');
rng(1);
[useboundary, displaydebug, useTemplateOptimizer ] = parseArguments(Porg,Qorg,K,pnorm,varargin);
radius = sqrt( ((max(Porg(1,:))-min(Porg(1,:)))^2  + (max(Porg(2,:))-min(Porg(2,:)) )^2)/4 );

%% Init constants
EPSILON=1e-3;
MMAX=50;
DMIN=0.00000001;
BPAD=0.15;
MAX_RETRIES=3;
initdelta=1;
%%
Nb=0;
N=size(Porg,2);
d=initdelta;
display(['initial delta: ' num2str(d)]);

% Boundary
if (useboundary)
    [Porg, Nb] = addBoundary(Porg,BPAD,N);
end
deformation.P=Porg;

P=Porg./radius;
Q=Qorg./radius;

% generate triangular mesh
TRI = delaunay(P(1,:),P(2,:));
ntri=size(TRI,1);
deformation.TRI=TRI;
% Prepare debug display
if (displaydebug)
    fig =figure ;
    energy=NaN+zeros(1,MMAX);
    coordinateChange=NaN+zeros(1,MMAX);
    Matchhistory = zeros(N,MMAX+1);
end

%Init frames at identity
Pdef = P;
W = (sum((Pdef(:,1:N)-Q).^2,1)'+d).^(pnorm/2 - 1);
W = W./max(W);
%     W=ones(N,1);                %all points have same weigth initially
Rotations=repmat([1 0 0 1]',1,ntri);% representing rotations with a,b: R=[a b;-b a]

myOptimizer = createOptimizationFunction(useTemplateOptimizer,N,Nb,P,Q,TRI,K);

m=0;
while (m<MMAX && d>DMIN)
    first=true;
    iter_per_delta=0;
    while (m<MMAX && (first || ~isfinite( beforeDeformationEnergy ) ...
            || (abs(beforeDeformationEnergy-afterDeformationEnergy) > EPSILON) ))
        
        %% convergence and sanity
        if (first)
            beforeDeformationEnergy = Inf;
        else
            beforeDeformationEnergy = afterDeformationEnergy;
        end
        
        if (exist( 'afterDeformationEnergy','var' ) && beforeDeformationEnergy-afterDeformationEnergy < -EPSILON)
            warning(['energy increased! ' num2str(beforeDeformationEnergy-afterDeformationEnergy)]);
            break;
        end
        
        %%
        m=m+1;
        iter_per_delta=iter_per_delta+1;
        
        %% Optimize
        oldPdef=Pdef;
        [Pdef, errorcode] = myOptimizer(Rotations,W);
        if errorcode ~= 0
            warning('Hmm, something went wrong with solver!');
            yalmiperror(errorcode)
            
            if errorcode == -1
                if ~exist('retry','var')
                    retry=1;
                end
                if retry<=MAX_RETRIES
                    display(['Retrying iteration... retry #' num2str(retry)]);
                    m=m-1;
                    Pdef=oldPdef;
                else
                    clear('retry');
                    error(['Retries failed on iteration #' m]);
                end
            end
        end
        
        %% Prepare for next iteration
        currCoordinateChange=norm(W'.*sum((Pdef(:,1:N)-oldPdef(:,1:N))),1)./(N*initdelta);
        afterDeformationEnergy=sum((sum((Pdef(:,1:N)-Q).^2,1)+d).^(pnorm/2))./N;
        
        W = (sum((Pdef(:,1:N)-Q).^2,1)'+d).^(pnorm/2 - 1);
        W = W./max(W);
        
        [Rotations, distortions, flips]=AnalyseDeformation(Pdef,TRI,P);
        
        if (displaydebug)
            Matchhistory(:,m+1)=W>0.5;
            energy(m)=afterDeformationEnergy;
            coordinateChange(m)=currCoordinateChange;
            displayDebug();
        end
        
        first=false;
    end
    
    d=d/2;
    display(['Updated delta: ' num2str(d)]);
end

%% Prepare Results
matches = find(W>0.5);
Pf=Porg(1:2,matches);
Qf=Qorg(1:2,matches);
deformation.Pdef=Pdef.*radius;
if (displaydebug)
    debugFigure=fig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END BDMatch
%% BEGIN inner functions
    function displayDebug()
        figure(fig);
        
        %         subplot(7,1,1);hist(W);
        %         subplot(7,1,2);plot(1:MMAX,energy);
        %         subplot(7,1,3);plot(1:MMAX,coordinateChange);
        %         subplot(7,1,4:5);imagesc(Matchhistory);
        %         subplot(7,1,6:7);
        
        subplot(2,4,1);hist(W); title('Weights');
        subplot(2,4,2);plot(1:MMAX,energy);title('Energy');
        subplot(2,4,6);plot(1:MMAX,coordinateChange); title('Coordinate change');
        subplot(2,4,5);imagesc(Matchhistory); title('Matching');
        subplot(2,4,[3 4 7 8]);
        
        hold off;
        %target
        scatter(Q(1,:),Q(2,:),25,[1 0 0],'filled');
        hold on
        %source
        scatter(Pdef(1,N+1:end),Pdef(2,N+1:end),8,[0 1 0],'filled');
        scatter(Pdef(1,1:N),Pdef(2,1:N),20,[0 0 1]);
        scatter(Pdef(1,W>0.9),Pdef(2,W>0.9),26,[0 1 1],'filled');
        %new mesh
        H=triplot(TRI,Pdef(1,:),Pdef(2,:));
        set(H(flips==1),'linestyle','--');
        set(gca,'ydir','reverse');
        axis equal;
        hold off;
        drawnow;
    end
%% END Inner functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function myOptimizer = createOptimizationFunction(useTemplateOptimizer,N,Nb,P,Q,TRI,K)
Pdefhat=sdpvar(2,N+Nb);
residuals=Pdefhat(:,1:N)-Q;
aux_r = sdpvar(2,length(residuals));
boundaryConstraints=[];

if (Nb>0)
    boundryA=sdpvar(2,2,'full');
    boundryt=sdpvar(2,1);
    
    boundaryConstraints= [Pdefhat(:,N+1:end) == boundryA*P(:,N+1:end) + repmat(boundryt,1,Nb)];
    
end

baseConstraints=[aux_r==residuals]+boundaryConstraints;

if (useTemplateOptimizer)
    templateOptimizer = createTemplateOptimizer(P, Pdefhat, residuals,...
        baseConstraints,TRI,K);
    myOptimizer = @(Rotations,W) (templateOptimizer{{Rotations',W'.^0.5}});
else
    solverOptions = sdpsettings('verbose',0,'cachesolvers',1,'solver','mosek');
    [faces] = GenerateFaces(P,TRI,Pdefhat);
    myOptimizer = @(Rotations,W) ...
        (sdpsolveOptimzer(Rotations,W, K, ...             % Optimization Parameters
        Pdefhat, residuals, baseConstraints, ... % Optimization Variables
        faces, solverOptions));         % helpers
end


end

function [faces] = GenerateFaces(P,TRI,Pdefhat)
tmp=sdpvar(2,2);
tface=struct('V',zeros(1,3),'A',tmp,'B',tmp(:),'C',tmp(:));

faces = repmat(tface,[size(TRI,1) 1]);

averager=eye(3)-(1/3)*ones(3,3);
for i=1:numel(faces)
    V=TRI(i,:);
    %             At =  Pdefhat(:,V)/([P(:,V); 1 1 1]);
    %             A=At(1:2,1:2);
    A=(Pdefhat(:,V)*averager)*pinv(P(:,V)*averager);
    B=(A-A'+trace(A)*eye(2))/2;
    C=(A+A'-trace(A)*eye(2))/2;
    faces(i)=struct('V', V, ...%                             'At',At,...
        'A', A,...
        'B',B(:), ...
        'C',C(:));
end

end

function [Pdef, errorcode] = sdpsolveOptimzer(Rotations,W, K, ...
    Pdefhat, residuals, baseConstraints, ...
    faces, solverOptions)
% create deformation constraints
ntri=length(faces);
N=size(residuals,2);
Constraints=baseConstraints;
%     for j=1:ntri
%         % rotate the space such that it includes the previous solution
%         % and more of the nonconvex space
%         R=[Rotations(2*j-1), Rotations(2*j); ...
%            -Rotations(2*j), Rotations(2*j-1)];
%         B = faces(j).B*R';
%         C = faces(j).C;%*Rotations{j}';
%
%         Constraints = Constraints + cone( C(:),((K-1)/(K+1))*(trace(B)/(sqrt(2))) );
%     end

Constraints = Constraints + generateBdConstraints(faces,Rotations,K);
% create objective via cone constraint.
t = sdpvar(1,1);
Constraints = Constraints + cone([sqrt(W').*residuals(1,:), sqrt(W').*residuals(2,:)], t);
Objective = 1000*t;

% solve
sol = solvesdp(Constraints,Objective,solverOptions);
if sol.problem == 0
    Pdef = double(Pdefhat);
else
    Pdef=zeros(size(Pdefhat));
end
errorcode = sol.problem;
end

%% parameterized optimizer .
% building a paramatrized optimizer object with 2 types of parameters, sqrtW and R.
function [optimizerTemplate, RotationParams, faces] = createTemplateOptimizer(P, Pdefhat, residuals, baseConstraints,TRI,K)%(N,Nb,P,Q,TRI,K)
N=size(residuals,2);
% Creating parametrized faces constraints
display('building constraints');
[faces Constraints RotationParams] = GenerateFacesConstraints(P,TRI,Pdefhat,K);
Constraints = Constraints + baseConstraints;

% build objective
aux_sqrtW=sdpvar(1,N); % Parameter for weights (sqrt(W_i))
aux_Obj=sdpvar(2,N);   % auxilary variable for residuals coordinates times sqrt(W_i)
% needed because parameters are not allowoed
% directly in cone
t=sdpvar(1,1);
display('building objective');

Constraints = Constraints + [aux_Obj == [aux_sqrtW;aux_sqrtW].*residuals];
Constraints = Constraints + cone(aux_Obj(:),t);
Objective = 1000*t;

% create template optimizer, forcing a mosek solver ('+')
options = sdpsettings('verbose',0,'cachesolvers',1,'solver','+mosek');
optimizerTemplate = optimizer(Constraints,Objective,options,{ RotationParams',aux_sqrtW},Pdefhat);
end

% generate parameterized constraints for all faces
function [faces, Constraints, RotationParams] = GenerateFacesConstraints(P,TRI,Pdefhat,K)
ntri=size(TRI,1);
faces = cell(ntri,1);

% 2 rotation parameters for ith face representing the Rotation
% part of the transformation. only 2 paramers per rotation are
% needed
RotationParams = sdpvar(2*ntri,1,'full');

% We use auxilary variables for the trace(B*R') because it includes
% variables multiplied by parameters which aren't allowed inside
% cone()
aux_trcBRt = sdpvar(ntri,1,'full');

averager=eye(3)-(1/3)*ones(3,3);

Constraints=[];
for i=1:ntri
    V=TRI(i,:);
    %         At =  Pdefhat(:,V)/([P(:,V); 1 1 1]);
    
    %         A=At(1:2,1:2);
    A=(Pdefhat(:,V)*averager)*pinv(P(:,V)*averager);
    B=(A-A'+trace(A)*eye(2))/2;
    C=(A+A'-trace(A)*eye(2))/2;
    
    %             faces{i}=struct('V', V, ...
    %                             'At',At,...
    %                             'A', A,...
    %                             'B',B, ...
    %                             'C',C);
    %
    
    Constraints = Constraints + [aux_trcBRt(i) == (RotationParams(2*i-1)*(B(1,1) + B(2,2)) - RotationParams(2*i)*(B(2,1)-B(1,2)))] ;
    Constraints = Constraints + [ cone(C(:), ((K-1)/(K+1))*(aux_trcBRt(i)/(sqrt(2)))) ];
end

end
%%

function [Rotations distortion flips ]=AnalyseDeformation(Pdef,TRI,P)
ntri =size(TRI,1);
Rotations = zeros(4,ntri);
distortion=zeros(ntri,1);
flips=zeros(ntri,1);
averager=eye(3)-(1/3)*ones(3,3);
for l=1:ntri
    V=TRI(l,:);
    %         At_ = Pdef(:,TRI(l,:))/([P(:,TRI(l,:)); 1 1 1]);
    %         A_=At_(1:2,1:2);
    A_=(Pdef(:,V)*averager)*pinv(P(:,V)*averager);
    [U, S, V, flip]=closest_rotation(A_);
    
    R=U*V';
    Rotations(:,l) = R(:);
    distortion(l)=abs(max(S(1,1),S(2,2))/min(S(1,1),S(2,2)));
    flips(l)=flip;
    
end

end
function [BdConstraints] = generateBdConstraints(faces,R,K)
ntri =length(faces);
BdConstraints=[];

% patch!
%             P=[Image.KeyPoints Image.BoundaryPoints];
%             tinyFaces = findTinyFaces(Image.TRI,P,0.00001);
%             faces=Image.faces(~tinyFaces);
%             R=R(:,~tinyFaces);
%end

C=[faces.C];
B=[faces.B];

BdConstraints = BdConstraints + ...
    cone([(B(1,:).*R(1,:) + B(3,:).*R(3,:) + B(2,:).*R(2,:) + B(4,:).*R(4,:)) * ... /trace(B*R')
    ((K-1)/((K+1)*sqrt(2))); C]) ;

end
function [P Nb] = addBoundary(P,PAD,N)
bb=[min(P(1,:)) max(P(1,:)) min(P(2,:)) max(P(2,:))];
boundry=boundrycoords(bb,PAD,N)';
Nb=size(boundry,2);
P=[P,boundry];
end

function [useboundary, displaydebug, useTemplateOptimizer] = parseArguments(P,Q,K,pnorm,varargin)
p=inputParser;
p.FunctionName='BDMatch';
p.addRequired('P',@(x) size(x,1)==2 && size(x,2)>1);
p.addRequired('Q',@(x) size(x,1)==2 && size(x,2)>1);
p.addRequired('K',@isscalar);
p.addRequired('pnorm',@isscalar);
p.addOptional('useboundary',false, @islogical)
p.addOptional('displaydebug',true, @islogical)
p.addOptional('useTemplateOptimizer',true, @islogical);

p.parse(P, Q, K , pnorm, varargin{1}{:});
pResults = p.Results;
useboundary = pResults.useboundary;
displaydebug = pResults.displaydebug;
useTemplateOptimizer = pResults.useTemplateOptimizer;
end