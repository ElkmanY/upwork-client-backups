clear all; clc
%% Pretreatment
Moon = imread('moon.png');
imshow(Moon)

MoonFil = stdfilt(Moon);
MoonBW = MoonFil > 1.2;
montage({Moon,MoonFil,MoonBW},'BackgroundColor','white');
%{
    Traversability map is obtained in MoonBW.
    0 (Black) : Crossable
    1 (White) : Uncrossable (obstacle)
%}

%% Dijkstra Algorithm
%Inits
%----------------------
sizeI = size(MoonBW);

% List of all crossable positions
Q = find(MoonBW==0);    % index array of all crossable positions
sizeQ = size(Q);

% capture of START & END position :

imshow(MoonBW)
hold on
title('A* shortest path | Arrival cost = ? ', 'FontSize', 20);
[x,y] = getpts ;    % get points by mouse
x = round(x); y = round(y);  % Round to integer value for indexing
x = min(max(x,1),sizeI(2)); y = min(max(y,1),sizeI(1)); % Limit to stay in image range
Pstart = [y(1), x(1)]; Pend = [y(2), x(2)];
hold off
close

%%%%% Uncomment this part to overide mouse capture and use a specific value for Pend %%%%%%%%%
% Pstart = [460 30]; Pend = [310 440]; % Can be modified without running disjkstra algorithm %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if selected positions are on crossable map
    % Index of Pstart & Pend in Image
    idxStart = sub2ind(sizeI, Pstart(1), Pstart(2));
    idxEnd = sub2ind(sizeI, Pend(1), Pend(2));

    % Checking positions are in crossable map
    if MoonBW(idxStart) == 1
        fprintf(2,'Start position is on an obstacle, Choose an other position ... \n')
        return
    end
    if MoonBW(idxEnd) == 1
        fprintf(2,'End position is on an obstacle, Choose an other position ... \n')
        return
    end
    
    
% Index of Pstart and Pend in Q
idxStartQ = find(Q == idxStart); 
idxEndQ = find(Q == idxEnd); 

% Initializing weights (All equals to Inf except Pstart = 0)
Weights = Inf * ones(sizeQ);
Weights(idxStartQ) = 0;

%Initialize Previous array with indefined elements
Previous = NaN * ones(sizeQ);

% ----------------------------
% Main Dijkstra Algorithm
% ----------------------------
tic % Start time counter

% Defining Neighbors Ofset and their weights (distance)
       % [(-1,-1)=sqrt(2)    (-1, 0)=1    (-1, 1=sqrt(2))
       %  ( 0,-1)=1             %P%       ( 0, 1)=1
       %  ( 1,-1)=sqrt(2)    ( 1, 0)=1    ( 1, 1)=sqrt(2)]
NeighborsOfs = [-1 -1; 0 -1; 1 -1; -1 0; 1 0; -1 1; 0 1; 1 1];
NeighborsWgt = [sqrt(2) 1 sqrt(2) 1 1 sqrt(2) 1 sqrt(2)]'; % Pythagorean theorem

Qtemp = Q;
WeightsTemp = Weights;
idxP = idxStart;
i = 0; % Show progress of the loop 
disp 'Dijkstra algorithm started ...';
while idxP ~= idxEnd%~isempty(Qtemp)
    i = i +1;
    if mod(i,1000) == 0
        i/1000
    end
    
    [MinWgt,idxMinWgt] = min(WeightsTemp);   % Value & Index of the minimum weight remaining
    idxP = Qtemp(idxMinWgt,:);          % Index of minimum weight in image
    [r,c] = ind2sub(sizeI,idxP);    % index to xy conversion
    P = [r,c];               % P takes xy coords of minimum weight

    
    % Tuning neighbors
    %
    Neighbors = P + NeighborsOfs;   % xy coords of the 8 neibors
    NeighborsLim = Neighbors((1<=Neighbors(:,1) & Neighbors(:,1)<=sizeI(1)) & ...
        (1<=Neighbors(:,2) & Neighbors(:,2)<=sizeI(2)),:); % To stay in image range
    NeighborsIdx = sub2ind(sizeI, NeighborsLim(:,1), NeighborsLim(:,2));    % xy to Idx Conversion
    
    [NeighborsQtemp,NeighborsIdxQt,ib] = intersect(Qtemp, NeighborsIdx); % Keep only neighbors existing in remaining Q 
    [~,NeighborsIdxQ,~] = intersect(Q, NeighborsQtemp); % take indices in Q of remaining neighbors
    NeighborsWgtQ = NeighborsWgt(ib); % Keep only weight ofsets of remaining Neighbors
    
    NeighborsNewDist = MinWgt+ NeighborsWgtQ; % Calculate new distance to the neighbors
    NeighborsNewDistIdx = NeighborsNewDist < Weights(NeighborsIdxQ); % idx of neighbors with min weight to be changed
    NeighborsNewDistIdxQ = NeighborsIdxQ(NeighborsNewDistIdx); % idx in Q of neighbors with min weight to be changed
    NeighborsNewDistIdxQt = NeighborsIdxQt(NeighborsNewDistIdx); % idx in Q of neighbors with min weight to be changed
    
    % Updating Weights and Previous array
    % Take new point as previous for all updated neighbors
    Previous(NeighborsNewDistIdxQ) = idxP; 
    % Update weights with minimum value
    Weights(NeighborsNewDistIdxQ) = NeighborsNewDist(NeighborsNewDistIdx);
    WeightsTemp(NeighborsNewDistIdxQt) = NeighborsNewDist(NeighborsNewDistIdx);    
    
    % Altering processed position (Positon and its weight)
    WeightsTemp(Qtemp == idxP) = [];    % Altering weight
    Qtemp(Qtemp == idxP) = [];          % Altering position
    

end
toc % End time counter and show time


% --------------------
% Drawing path
% --------------------

% Based on 'Previous', we start a reverse loop from END to START to extract
% the shortest path

ShortPath = idxEnd; % Array to fill with shortest path
j = idxEnd;

while j ~= idxStart
    idxjQ = find(Q == j); 
    j = Previous(idxjQ);
    ShortPath = [j, ShortPath];
end

[PathRow, PathCol] = ind2sub(sizeI,ShortPath); %index to xy coords conversion
PathXY = [PathCol' PathRow'];

Distance = Weights(Q == idxEnd); %Distance from START to END is END weight
save('Path_Dijkstra.mat','Previous','ShortPath','PathXY','Weights')

% Drawing shortest past in filtred image
imshow(MoonFil);
hold on;
scatter(Pstart(2),Pstart(1),'g','LineWidth',5); % Draw green point for START
scatter(Pend(2),Pend(1),'b','LineWidth',5);     % Draw blue point  for END
plot(PathCol, PathRow, 'r*', 'LineWidth', 0.5, 'MarkerSize', 2); % Draw path on red
title(['Dijkstra shortest path | Arrival cost = ', num2str(Distance)], 'FontSize', 20);
hold off;
