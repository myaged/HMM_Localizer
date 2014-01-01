clear, close all

%% HMM Robot localization
% M. Yagci Date: 23.12.2013

% -------------------------------------------------------------------
%% Map initialization + Observations
% -------------------------------------------------------------------
nRows = 16; nCols = 4; 

obstacles = [ 	1 2; 1 3; 2 3; 3 1; 5 2; 5 3; 5 4; 7 1; 7 2; 7 3; 8 2; 8 3; 
                10 3; 11 4; 12 1; 12 3; 14 2; 14 3; 15 2; 15 3; 15 4; 16 3];
sizeObstacles = size(obstacles,1);

% state neighborhood model
nDirs = 4; % (N,W,E,S) 
neighborhood = zeros(nRows,nCols,nDirs);

pdf = zeros(nRows,nCols);
for i = 1:nRows
    for j=1:nCols
		pdf(i,j) = 1/(nRows*nCols - sizeObstacles);
		
		% mark obstacles for each neighborhood direction, e.g. [1 1 0 0] 
		if ( j+1 > nCols || sum(ismember(obstacles,[i j+1],'rows'))==1)
			neighborhood(i,j,1) = 1;
		end
		if ( j-1 < 1 || sum(ismember(obstacles,[i j-1],'rows'))==1)
			neighborhood(i,j,4) = 1;
		end
		if ( i+1 > nRows || sum(ismember(obstacles,[i+1 j],'rows'))==1)
			neighborhood(i,j,3) = 1;
		end
		if ( i-1 < 1 || sum(ismember(obstacles,[i-1 j],'rows'))==1)
			neighborhood(i,j,2) = 1;		
		end		
    end
end

% sensor error
err = 0.2; 

% observed obstacles in directions N,W,E,S
observations = [1 1 0 1; 
 				1 0 0 1;
 				1 0 0 0;
 				0 1 0 0];
%               0 0 1 0
%                 ];
				

% -------------------------------------------------------------------
%% Forward algorithm for filtering
% -------------------------------------------------------------------

% Initialize forward variables
for i = 1:nRows
	for j=1:nCols	
		if (sum(ismember(obstacles,[i j],'rows')) == 0)
			nv = [ neighborhood(i,j,1) neighborhood(i,j,2) neighborhood(i,j,3) neighborhood(i,j,4) ];
			observations(1,:);
			d = sum(abs(nv-observations(1,:)));
			pdf(i,j) = pdf(i,j) * (1-err)^(4-d)*err^d;
		else
			pdf(i,j) = 0;
		end
	end	
end
alphas = pdf;

% Recursion
for o=2:size(observations,1)

	currentObservation = observations(o,:);

	tempGrid = zeros(nRows,nCols); % holds intermediate scores to calculate forward variables (alphas)
	for i = 1:nRows
		for j=1:nCols
			if (sum(ismember(obstacles,[i j],'rows')) == 0)
				nv = [ neighborhood(i,j,1) neighborhood(i,j,2) neighborhood(i,j,3) neighborhood(i,j,4) ];

				pTrans = sum(nv)/length(nv);

				if nv(1) == 0
					tempGrid(i,j+1) = ...
                        tempGrid(i,j+1) + pTrans * alphas(i,j);
				end
				if nv(2) == 0
					tempGrid(i-1,j) = ...
                        tempGrid(i-1,j) + pTrans * alphas(i,j);
				end
				if nv(3) == 0
					tempGrid(i+1,j) = ...
                        tempGrid(i+1,j) + pTrans * alphas(i,j);
				end
				if nv(4) == 0
					tempGrid(i,j-1) = ...
                        tempGrid(i,j-1) + pTrans * alphas(i,j);
				end
			end			
		end
	end
	
	for i = 1:nRows
		for j=1:nCols
			if (sum(ismember(obstacles,[i j],'rows')) == 0)
				nv = [ neighborhood(i,j,1) neighborhood(i,j,2) neighborhood(i,j,3) neighborhood(i,j,4) ];
				d = sum(abs(nv-currentObservation));
				b = (1-err)^(4-d)*err^d;
				alphas(i,j) = tempGrid(i,j) * b;
			else
				alphas(i,j) = 0;
			end
		end
	end
end

% normalized alphas to get the posterior
posterior = alphas/ sum(sum(alphas));


% -------------------------------------------------------------------
%% Plot map with posterior
% -------------------------------------------------------------------
sizeGridX = nRows+1;
sizeGridY = nCols+1;
for i=1:sizeGridY
    plot([1 sizeGridX]-0.5,[i i]-0.5,'-','color',[.6,.6,.6]);hold on;
end
for i=1:sizeGridX
    plot([i i]-0.5,[1 sizeGridY]-0.5,'-','color',[.6,.6,.6]);hold on;
end
axis off

for i = 1:nRows
    for j=1:nCols
        isObstacle = false;
        for o=1:sizeObstacles
            if i == obstacles(o,1) && j == obstacles(o,2)
                isObstacle = true; break;
            end
        end
        
        if isObstacle
            text(i-.28,j,'X','color',[.6 .6 .6], ...
                'fontsize', 20, 'BackgroundColor',[.6 .6 .6])

                
        else
        	ms = posterior(i,j)*20 +.00000001;
        	plot(i,j,'marker', 'o', 'markerfacecolor', 'black', 'markersize', ms)
        end
    end
end

hFig = figure(1);
set(hFig, 'Position', [100 100 700 250])
