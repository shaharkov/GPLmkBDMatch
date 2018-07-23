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


function [ matches, W, deformation, debugFigure ] = BDMatchMODIFIED(leftF, leftV, leftInd, leftBoundaryInd, rightCoords, K ,pnorm, initdelta, boundaryConstraintType, displaydebug)

%% Init + constants
if ~exist('displaydebug','var')
    displaydebug = true;
end
EPSILON=1e-3;
MMAX=50;
DMIN=0.00000001;
MAX_RETRIES=3;


%%
d=initdelta;
display(['initial delta: ' num2str(d)]);

N = numel(leftInd);

% Prepare debug display
if (displaydebug)
    energy=nan(1,MMAX);
    coordinateChange=nan(1,MMAX);
    Matchhistory = zeros(N,MMAX+1);
end

%Init frames at identity
ntri = size(leftF,2);
leftVdef = leftV;
residuals=leftVdef(:,leftInd)-rightCoords;
W = (sum(residuals.^2,1)'+d).^(pnorm/2 - 1);
W = W./max(W);
%     W=ones(N,1);                %all points have same weigth initially
[Rotations, distortions, flips]=AnalyseDeformation(leftVdef,leftF',leftV);

if (displaydebug)
    fig =figure;
    displayDebug();
    displayDebug();
    fig =figure;
end

myOptimizer = createOptimizationFunction(leftF, leftV, leftInd, leftBoundaryInd, rightCoords, K, boundaryConstraintType);

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
        oldleftVdef=leftVdef;
        [leftVdef, errorcode] = myOptimizer(Rotations,W);
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
                    leftVdef=oldleftVdef;
                else
                    clear('retry');
                    error(['Retries failed on iteration #' m]);
                end
            end
        end
        
        %% Prepare for next iteration
        currCoordinateChange=norm(W'.*sum((leftVdef(:,leftInd)-oldleftVdef(:,leftInd))),1)./(N*initdelta);
        residuals=leftVdef(:,leftInd)-rightCoords;
        afterDeformationEnergy=sum((sum(residuals.^2,1)+d).^(pnorm/2))./N;
        
        W = (sum(residuals.^2,1)'+d).^(pnorm/2 - 1);
        W = W./max(W);
        
        [Rotations, distortions, flips]=AnalyseDeformation(leftVdef,leftF',leftV);
        
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
deformation.leftVdef = leftVdef;
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
        
        %hold off;
        cla
        patch('vertices',leftVdef','faces',leftF','facecolor','none','edgecolor',0.9*[1 1 1])
        hold on
        line([rightCoords(1,:); leftVdef(1,leftInd)], [rightCoords(2,:); leftVdef(2,leftInd)],'color',[0.8 1 0.8]);
        line([rightCoords(1,(W>0.9)); leftVdef(1,leftInd(W>0.9))], [rightCoords(2,(W>0.9)); leftVdef(2,leftInd(W>0.9))],'color',[1 0.8 0.8]);

%         %target
%         scatter(rightCoords(1,:),rightCoords(2,:),25,[1 0 0],'filled');
%         %source
%         scatter(leftVdef(1,leftBoundaryInd),leftVdef(2,leftBoundaryInd),8,[0.8 0 0.8],'filled');
%         
%         scatter(leftVdef(1,leftInd(W>0.9)),leftVdef(2,leftInd(W>0.9)),26,[0 0.5 1],'filled');
%         scatter(leftVdef(1,leftInd),leftVdef(2,leftInd),20,[0 0 1]);

        %target
        scatter(rightCoords(1,:),rightCoords(2,:),25,[1 0 0]);
        %source
        scatter(leftVdef(1,leftBoundaryInd),leftVdef(2,leftBoundaryInd),8,[0.8 0 0.8],'filled');
        
        scatter(leftVdef(1,leftInd(W>0.9)),leftVdef(2,leftInd(W>0.9)),26,[1 0 0],'filled');
        scatter(leftVdef(1,leftInd),leftVdef(2,leftInd),20,[0 0 1]);
        
        %new mesh
        %         H=triplot(leftF',leftVdef(1,:),leftVdef(2,:));
        %         set(H(flips==1),'linestyle','--');
        
        
        %set(gca,'ydir','reverse');
        axis equal;
        %hold off;
        drawnow;
    end
%% END Inner functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function myOptimizer = createOptimizationFunction(leftF, leftV, leftInd, leftBoundaryInd, rightCoords, K, boundaryConstraintType)
leftVdefhat=sdpvar(2, size(leftV,2));
residuals=leftVdefhat(:,leftInd)-rightCoords;

switch boundaryConstraintType
    case 'none'
        boundaryConstraints = [];
    case 'fixed'
        boundaryConstraints = [leftVdefhat(:,leftBoundaryInd)==leftV(:,leftBoundaryInd)];
    case 'linear'
        boundryA=sdpvar(2,2,'full');
        boundaryConstraints= [leftVdefhat(:,leftBoundaryInd)==boundryA * leftV(:,leftBoundaryInd)];
    case 'affine'
        boundryA=sdpvar(2,2,'full');
        boundryt=sdpvar(2,1);
        boundaryConstraints= [leftVdefhat(:,leftBoundaryInd)==boundryA * leftV(:,leftBoundaryInd) + repmat(boundryt,[1 length(leftBoundaryInd)])];
    case 'similarity'
        a = sdpvar(2,1);
        boundryA=[a(1), a(2); -a(2), a(1)];
        boundaryConstraints= [leftVdefhat(:,leftBoundaryInd)==boundryA * leftV(:,leftBoundaryInd)];
    case 'similarity+translation'
        a = sdpvar(2,1);
        boundryt=sdpvar(2,1);
        boundryA=[a(1), a(2); -a(2), a(1)];
        boundaryConstraints= [leftVdefhat(:,leftBoundaryInd)==boundryA * leftV(:,leftBoundaryInd) + repmat(boundryt,[1 length(leftBoundaryInd)])];
        
    otherwise
        error('invalid boundaryConstraintType');
end



% % if (Nb>0)
% %     boundryA=sdpvar(2,2,'full');
% %     boundryt=sdpvar(2,1);
% %     
% %     boundaryConstraints= [Pdefhat(:,N+1:end) == boundryA*P(:,N+1:end) + repmat(boundryt,1,Nb)];
% % end

baseConstraints=boundaryConstraints;

solverOptions = sdpsettings('verbose',0,'cachesolvers',1,'solver','mosek');
[faces] = GenerateFaces(leftV,leftF',leftVdefhat);
myOptimizer = @(Rotations,W) ...
    (sdpsolveOptimzer(Rotations,W, K, ...             % Optimization Parameters
    leftVdefhat, residuals, baseConstraints, ... % Optimization Variables
    faces, solverOptions));         % helpers
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

function [leftVdef, errorcode] = sdpsolveOptimzer(Rotations,W, K, ...
    leftVdefhat, residuals, baseConstraints, ...
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
    leftVdef = double(leftVdefhat);
else
    leftVdef=zeros(size(leftVdefhat));
end
errorcode = sol.problem;
end


%%
function [Rotations, distortion, flips ]=AnalyseDeformation(leftVdef,leftF,leftV)
ntri =size(leftF,1);
Rotations = zeros(4,ntri);
distortion=zeros(ntri,1);
flips=zeros(ntri,1);
averager=eye(3)-(1/3)*ones(3,3);
for l=1:ntri
    V=leftF(l,:);
    %         At_ = Pdef(:,TRI(l,:))/([P(:,TRI(l,:)); 1 1 1]);
    %         A_=At_(1:2,1:2);
    A_=(leftVdef(:,V)*averager)*pinv(leftV(:,V)*averager);
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
    ((K-1)/((K+1)*sqrt(2))); C]);

BdConstraints = BdConstraints + ...
    cone([(B(1,:).*R(1,:) + B(3,:).*R(3,:) + B(2,:).*R(2,:) + B(4,:).*R(4,:)) * ... /trace(B*R')
    (1/sqrt(2))-0.75; C]); % restrict scale to be not less than 0.75

end