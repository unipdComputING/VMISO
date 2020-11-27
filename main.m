clear all; clear; clc;
%---------------------------------------------------------------------INPUT
E = 210000; %young modulus MPa
nu = 0.1; %Poisosn's ration

sy = TABLE(1);
%Plastic moduli
sy.data=([   0.000, 200.0;
             0.001, 220.0;
             0.010, 200.0;
             0.015, 180.0;
             0.1,   150.0;]);
% sy.data=([   0.000, 200.0;
%              0.001, 250.0;
%              0.010, 300.0;]);

dstrain = [-0.00005; -0.00005; 0.0001; 0.0; 0.0; 0.0;];
%dstrain = [0.0; 0.0; 0.002; 0.0; 0.0; 0.0;];
TOTINC = 210;
LOADUNLOAD = [30;30;30;30;30];
%--------------------------------------------------------------------------
strain = zeros(6, 1);

STATEV   = zeros(10, 1);
%note: the real statev dimension is 7 [stress; alpha] but this vector
%      has been used also to store some interesting variables as:
%      the Von-Mises equivalent stress and the equivalent strain
STRESS   = zeros(6, TOTINC);
STRAIN   = zeros(6, TOTINC);
PSTRAIN  = zeros(6, TOTINC);
GAMMA    = zeros(1, TOTINC);
SVM      = zeros(1, TOTINC);
EQSTRAIN = zeros(1, TOTINC);
mat = VM3D(E, nu, sy);
check = true;
lu = LOADUNLOAD(1);
lu_cont = 1;
lu_tot = size(LOADUNLOAD, 1);
for i=1:TOTINC
  fprintf('++INC: %i\n', i);
  [stress, D, STATEV] = mat.GETCONSTMAT(strain, dstrain, STATEV);
  strain = strain + dstrain;
  if ((i == lu) && (lu_cont < lu_tot))
    dstrain = -1 * dstrain;
    lu_cont = lu_cont + 1;
    lu = lu + LOADUNLOAD(lu_cont);
  end

  STRESS  (:, i) = stress(:);
  STRAIN  (:, i) = strain(:);
  PSTRAIN (:, i) = STATEV(1:6);
  GAMMA   (1, i) = STATEV(7);
  SVM     (1, i) = STATEV(8);
  EQSTRAIN(1, i) = STATEV(9);
  check          = STATEV(10);
  if (check == false)
    break
  end
end
%plot(EQSTRAIN, SVM,'k-');
plot(STRAIN(3,:), STRESS(3,:),'k-');
%--------------------------------------------------------------------------

