% Calculate homogenous joint transformation matrices for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% T_mdh [4x4x6]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(6+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 18:09
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = palh2m2OL_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 18:02:16
% EndTime: 2020-06-30 18:02:16
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->8), mult. (0->0), div. (0->0), fcn. (24->12), ass. (0->13)
t49 = cos(qJ(1));
t48 = cos(qJ(2));
t47 = cos(qJ(3));
t46 = cos(qJ(4));
t45 = cos(qJ(5));
t44 = cos(qJ(6));
t43 = sin(qJ(1));
t42 = sin(qJ(2));
t41 = sin(qJ(3));
t40 = sin(qJ(4));
t39 = sin(qJ(5));
t38 = sin(qJ(6));
t1 = [t49, -t43, 0, 0; t43, t49, 0, 0; 0, 0, 1, 0; t48, -t42, 0, pkin(1); 0, 0, -1, 0; t42, t48, 0, 0; t47, -t41, 0, pkin(4); t41, t47, 0, 0; 0, 0, 1, 0; t46, -t40, 0, pkin(2); t40, t46, 0, 0; 0, 0, 1, 0; t45, -t39, 0, pkin(5); t39, t45, 0, 0; 0, 0, 1, 0; t44, -t38, 0, pkin(3); 0, 0, 1, 0; -t38, -t44, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,6);             % numerisch
else,                         T_mdh = sym('xx', [4,4,6]); end % symbolisch

for i = 1:6
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
