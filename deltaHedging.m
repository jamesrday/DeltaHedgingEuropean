% The following program delta hedges for a given long european call option

clear all
clc

% Initialise variables

iT = 0.25; % Time to maturity in years
iSteps = 1000; % Number of steps (except initial step t = 0)
dt = iT./iSteps; % Time interval length
iHedgeFreq = dt; % Time interval for hedging ('dt' is quasi-continuous)
iNumSimu = 100; % Number of simulations
ir = 0.005; % Risk free rate
iS0 = 100; % Initial price of the stock
iK = 100; % Option strike
iMu = 0.05; % Expected return of the stock
iSigmaImp = 0.2; % Implied volatility of the stock
iSigmaAct = 0.3; % Actual volatility of the stock
iCost = 0; % Transaction costs

% Define paths of underlying

Time = [0:dt:iT]'; % Time vector
StockPaths = zeros(iSteps + 1, iNumSimu);
StockPaths(1,:) = iS0;
Sigma = iSigmaAct; 

for i = 1 : iNumSimu
    for j = 2 : iSteps+1
        dS = iMu*StockPaths(j-1,i)*dt + Sigma*StockPaths(j-1,i)*sqrt(dt)*randn();
        StockPaths(j,i) = StockPaths(j-1,i) + dS;
    end
end

% Calculation of option price and delta for each time step (observed)

OptPriceAct = zeros(iSteps + 1, iNumSimu);
OptDeltaAct = zeros(iSteps + 1, iNumSimu);

for i = 1 : iNumSimu
    for j = 1 : iSteps + 1
        d1 = (log(StockPaths(j,i)./iK) + (ir + 0.5 .* (iSigmaAct^2)) .* (iT - Time(j))) ./ (iSigmaAct * sqrt(iT - Time(j)));
        d2 = (log(StockPaths(j,i)./iK) + (ir - 0.5 .* (iSigmaAct^2)) .* (iT - Time(j))) ./ (iSigmaAct * sqrt(iT - Time(j)));
        
        OptDeltaAct(j,i) = normcdf(d1);
        OptPriceAct(j,i) = StockPaths(j,i) .* OptDeltaAct(j,i) - iK .* exp(-ir .* (iT - Time(j))) .* normcdf(d2);
    end
end

% Calculation of option price and delta for each time step (impled)

OptPriceImp = zeros(iSteps + 1, iNumSimu);
OptDeltaImp= zeros(iSteps + 1, iNumSimu);

for i = 1 : iNumSimu
    for j = 1 : iSteps
        d1 = (log(StockPaths(j,i)./iK) + (ir + 0.5 .* (iSigmaImp^2)) .* (iT - Time(j))) ./ (iSigmaImp * sqrt(iT - Time(j)));
        d2 = (log(StockPaths(j,i)./iK) + (ir - 0.5 .* (iSigmaImp^2)) .* (iT - Time(j))) ./ (iSigmaImp * sqrt(iT - Time(j)));
        
        OptDeltaImp(j,i) = normcdf(d1);
        OptPriceImp(j,i) = StockPaths(j,i) .* OptDeltaImp(j,i) - iK .* exp(-ir .* (iT - Time(j))) .* normcdf(d2);
    end
end

% Observed Volatility Hedging

oPortfolioAct = zeros(iSteps + 1, iNumSimu);
oCashAct = zeros(iSteps + 1, iNumSimu);
oPosiAct = zeros(iSteps + 1, iNumSimu); % total position
oPnLAct = zeros(iSteps + 1, iNumSimu);
oPortfolioAct(1,:) = OptPriceImp(1,:) - OptDeltaAct(1,:) .* StockPaths(1,:); % Long the cheap option and short the hedge (short the expensive "option")
oCashAct(1,:) = -oPortfolioAct(1,:); % Initial cash is the different between the cost of the option and the cash obtaind afeter short selling the stock
oPosiAct(1,:) = oPortfolioAct(1,:) + oCashAct(1,:); % Initial position
oPnLAct(1,:) = oPosiAct(1,:) - abs(OptDeltaAct(1,:)) .* iCost; % Initial profit/loss is zero + transaction costs
oNumHedges = 0;

for i = 1 : iNumSimu
    Delta = OptDeltaAct(1,i); % Initial delta remains constant until the first re-hedge takes place
    
    for j = 2 : iSteps+1
        if rem(Time(j), iHedgeFreq) == 0 
            oPortfolioAct(j,i) =  OptPriceImp(j,i) - OptDeltaAct(j,i) .* StockPaths(j,i);
            oCashAct(j,i) = oCashAct(j-1,i).*(1 + ir.*dt) + (OptDeltaAct(j,i) - Delta) .* StockPaths(j,i); % Subsequent cash positions are the sum of the preceeding cash postion + interest income(+/-) - cash needed/obtained to adjust the delta hedge
            oPosiAct(j,i) = oPortfolioAct(j,i) + oCashAct(j,i);
            oPnLAct(j,i) = oPnLAct(j-1,i) + oPosiAct(j,i) - oPosiAct(j-1,i) - abs(OptDeltaAct(j,i) - Delta) .* iCost;
            
            Delta = OptDeltaAct(j,i);
            oNumHedges = oNumHedges + 1;
        else
            oPortfolioAct(j,i) = OptPriceImp(j,i) - Delta.*StockPaths(j,i);
            oCashAct(j,i) = oCashAct(j-1,i).*(1 + ir.*dt);
            oPosiAct(j,i) = oPortfolioAct(j,i) + oCashAct(j,i);
            oPnLAct(j,i) = oPnLAct(j-1,i) + oPosiAct(j,i) - oPosiAct(j-1,i);
        end
    end
end

oPnLAct(iSteps + 1, :) = oPnLAct(iSteps + 1, :) + max(StockPaths(iSteps + 1,:) - iK, 0); % At expiration, if the option is OTM it is executed

oRealOptionAct = OptPriceImp;
oArtOptionAct = -1.*(oPortfolioAct - oRealOptionAct + oCashAct);

% Display observed profit and loss
plot(oPnLAct)

% Implied Volatility Hedging

%oPortfolioImp = zeros(iSteps + 1, iNumSimu);
%oCashImp = zeros(iSteps + 1, iNumSimu);
%oPosiImp = zeros(iSteps + 1, iNumSimu); % total position
%oPnLImp = zeros(iSteps + 1, iNumSimu);
%oPortfolioImp(1,:) = OptPriceImp(1,:) - OptDeltaImp(1,:) .* StockPaths(1,:); % Long the cheap option and short the hedge (short the expensive "option")
%oCashImp(1,:) = -oPortfolioImp(1,:); % Initial cash is the different between the cost of the option and the cash obtaind afeter short selling the stock
%oPosiImp(1,:) = oPortfolioImp(1,:) + oCashImp(1,:); % Initial position
%oPnLImp(1,:) = oPosiImp(1,:) - abs(OptDeltaImp(1,:)) .* iCost; % Initial profit&loss is zero + transaction costs
%oNumHedges = 0;

%for i = 1 : iNumSimu
%    Delta = OptDeltaImp(1,i); % Initial delta remains constant until the first re-hedge takes place
%    for j = 2 : iSteps+1
%        if rem(Time(j), iHedgeFreq) == 0
%            oPortfolioImp(j,i) =  OptPriceImp(j,i) - OptDeltaImp(j,i) .* StockPaths(j,i);
%            oCashImp(j,i) = oCashImp(j-1,i).*(1 + ir.*dt) + (OptDeltaImp(j,i) - Delta) .* StockPaths(j,i); % Subsequent cash positions are the sum of the preceeding cash postion + interest income(+/-) - cash needed/obtained to adjust the delta hedge
%            oPnLImp(j,i) = oPnLImp(j-1,i) + oPosiImp(j,i) - oPosiImp(j-1,i) - abs(OptDeltaImp(j,i) - Delta) .* iCost;
            
%            Delta = OptDeltaImp(j,i);
%            oNumHedges = oNumHedges + 1;
%        else
%            oPortfolioImp(j,i) = OptPriceImp(j,i) - Delta.*StockPaths(j,i);
%            oCashImp(j,i) = oCashImp(j-1,i).*(1 + ir.*dt);
%            oPnLImp(j,i) = oPnLImp(j-1,i) + oPosiImp(j,i) - oPosiImp(j-1,i);
%        end
%    end    
%end

%oPnLImp(iSteps + 1, :) = oPnLImp(iSteps + 1, :) + max(StockPaths(iSteps + 1,:) - iK, 0); % At expiration, if the option is OTM it is executed

%oRealOptionImp = OptPriceImp;
%oArtOptionImp = -1.*(oPortfolioImp - oRealOptionImp + oCashImp);

%plot(oPnLImp)

