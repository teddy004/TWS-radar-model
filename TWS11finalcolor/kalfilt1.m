function [filtered, residuals , covariances, kalmgains,xhat,range1,Azimuth1,Elevation1] = kalfilt1(trajectory, xhatminus, Pminus, phi, R, Q,Azi,El )
% close all;
% clear all;
% clc all;

% USAGE: [filtered, residuals , covariances, kalmgains] = kalfilt(trajectory, x0, P0, phi, R, Q )
%
% INPUTS
%
% name         dimension                    explanation                                  units
%------        ------                       ---------------                              -------
% trajectory   NUMMEASUREMENTS X NUMPOINTS  trajectory in radar reference coords         [m;m;m]
% x0           NUMSTATES X 1                initial estimate of state vector             m, m/s
% P0           NUMSTATES X NUMSTATES        initial estimate of covariance matrix        m, m/s
% phi          NUMSTATES X NUMSTATES        state transition matrix                      -
% R            NUMMEASUREMENTS X NUMMEASUREMENTS   measurement error covariance matrix   m
% Q            NUMSTATES X NUMSTATES        state error covariance matrix                m, m/s
%
% OUTPUTS
%
% name         dimension                    explanation                                  units
%------        ------                       ---------------                              -------
% filtered     NUMSTATES X NUMPOINTS        filtered trajectory x,y,z pos, vel    [m; m/s; m; m/s; m; m/s]
% residuals    NUMSTATES X NUMPOINTS        residuals of filtering                [m;m;m]
% covariances  NUMSTATES X NUMPOINTS        diagonal of covariance matrix         [m;m;m]
% kalmgains    (NUMSTATES X NUMMEASUREMENTS) 
%                 X NUMPOINTS               Kalman gain matrix                    -
%

%%
NUMSTATES = 6 ;
NUMMEASUREMENTS = 3 ;
NUMPOINTS = size(trajectory, 2) ;
%%
% initialize output matrices
filtered = zeros(NUMSTATES, NUMPOINTS) ;
xhatminus_next = zeros(NUMSTATES, NUMPOINTS) ;

residuals = zeros(NUMMEASUREMENTS, NUMPOINTS) ;
Y = zeros(NUMMEASUREMENTS, NUMPOINTS);
covariances = zeros(NUMSTATES, NUMPOINTS) ; 
kalmgains = zeros(NUMSTATES*NUMMEASUREMENTS, NUMPOINTS) ;
%%
% set matrix relating measurements to states
H = [1 0 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 0 1 0];
global pminus_n
global xp_n


   %%
for loop = 1: NUMPOINTS

   
   % compute the Kalman gain
   S = H*Pminus*H' + R;
   K = Pminus*H'*inv(S);
   kalmgains(:,loop) = reshape(K, NUMSTATES*NUMMEASUREMENTS, 1) ;
   
   % update the estimate with the measurement z
   z = trajectory(:,loop);
   m=sqrt(z(1)^2 + z(2)^2 + z(3)^2);
   residuals(:,loop)=z - H*xhatminus;
   Y = residuals(:,loop);
   d=sqrt(Y'*inv(S)*Y);
   w=H*xhatminus;
   range=sqrt(w(1)^2 + w(2)^2 + w(3)^2);
   [costMat] = gattting1(range,d,m);
   [R] = munkres1(costMat);
   xt=real(R*sin(El)*cos(Azi));
   yt=real(R*sin(El)*sin(Azi));
   zt=real(R*sin(El));
   xhatminus=[xt,300,yt,300,zt,300]';
   
   
   %UPDATE THE VALUE
   xhat = xhatminus + K*(z - H*xhatminus);
   position=H*xhat;
   range1=real(sqrt(position(1)^2 + position(2)^2 + position(3)^2));
   Azimuth1=atan(position(2)/position(1));
        
   Elevation1=atan(position(3)/(sqrt(position(1)^2 + position(2)^2)));
   filtered(:,loop) = xhat ;
   residuals(:,loop)=real(z - H*xhatminus);
   %residuals(:,loop) = xhat - xhatminus ;
   % update the error covariance for the updated estimate
   P =real( ( eye(NUMSTATES, NUMSTATES) - K*H)*Pminus) ;
   covariances(:,loop) = diag(P) ;  % only save diagonal of covariance matrix
   
   % project ahead
   xhatminus_next = phi*xhat ;
   
   Pminus_next = phi*P*phi' + Q ;
   
   xhatminus = xhatminus_next ;
   pminus_n= Pminus_next ;
   xp_n = xhatminus_next;
   
end % next loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
