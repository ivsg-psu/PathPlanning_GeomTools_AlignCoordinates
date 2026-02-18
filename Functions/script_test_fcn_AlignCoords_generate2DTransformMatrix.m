% script_test_fcn_AlignCoords_generate2DTransformMatrix
% Tests function: fcn_AlignCoords_generate2DTransformMatrix

% 
% REVISION HISTORY:
% 
% 2023_03_23 by Sean Brennan
% - first write of function


close all;



S = 2;
theta = 70*pi/180;
tx = 2;
ty = 7;
Tscale = [...
         S                   0                 0; 
         0                   S                 0; 
         0                   0                 1 
    ];
Trotate = [...
     cos(theta)  -sin(theta) 0; 
      sin(theta)  cos(theta) 0; 
         0                   0                 1 
    ];

Ttranslate = [...
         1                   0                 tx; 
         0                   1                 ty; 
         0                   0                 1 
    ];

%% Try basic call
figNum = 1;
order_string = 'rts';
T = fcn_AlignCoords_generate2DTransformMatrix( ...
    S, theta, tx, ty, order_string, figNum);

Trts = Trotate*Ttranslate*Tscale;
Ttrue = Trts;
assert(isequal(T,Ttrue));

%% Try another
figNum = 2;
order_string = 'srt';
T = fcn_AlignCoords_generate2DTransformMatrix( ...
    S, theta, tx, ty, order_string, figNum);
Tsrt = Tscale*Trotate*Ttranslate;
Ttrue = Tsrt;
assert(isequal(T,Ttrue));

%% Fail conditions
if 1==0
    %% Bad input - wrong number of arguments
    fcn_AlignCoords_generate2DTransformMatrix(24);

    %% Bad input - wrong number of arguments
    fcn_CodeX_generateNumbersLike_KEEP(24,3,4);
    %% Bad input - not integer
    fcn_CodeX_generateNumbersLike_KEEP(2.34,5);
    %% Bad input - not positive
    fcn_CodeX_generateNumbersLike_KEEP(-2,5);
    %% Bad input - not integer
    fcn_CodeX_generateNumbersLike_KEEP(2,5.4);
    %% Bad input - not positive
    fcn_CodeX_generateNumbersLike_KEEP(2,-5.4)
end
