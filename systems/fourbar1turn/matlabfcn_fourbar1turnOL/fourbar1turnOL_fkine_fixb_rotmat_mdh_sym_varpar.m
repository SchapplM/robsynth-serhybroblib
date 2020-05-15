% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = fourbar1turnOL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:44
% EndTime: 2020-04-12 19:40:44
% DurationCPUTime: 0.10s
% Computational Cost: add. (63->36), mult. (40->30), div. (0->0), fcn. (84->10), ass. (0->23)
t13 = sin(qJ(4));
t15 = sin(qJ(1));
t23 = t15 * t13;
t17 = cos(qJ(2));
t22 = t15 * t17;
t18 = cos(qJ(1));
t21 = t18 * t13;
t20 = t18 * t17;
t12 = qJ(2) + qJ(3);
t11 = pkin(5) + 0;
t14 = sin(qJ(2));
t19 = t14 * pkin(2) + t11;
t16 = cos(qJ(4));
t9 = qJ(5) + t12;
t8 = cos(t12);
t7 = sin(t12);
t6 = t16 * pkin(4) + pkin(1);
t5 = cos(t9);
t4 = sin(t9);
t3 = t18 * t16;
t2 = t15 * t16;
t1 = t17 * pkin(2) - pkin(3) * t8;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t18, -t15, 0, 0; t15, t18, 0, 0; 0, 0, 1, t11; 0, 0, 0, 1; t20, -t18 * t14, t15, 0; t22, -t15 * t14, -t18, 0; t14, t17, 0, t11; 0, 0, 0, 1; -t18 * t8, t18 * t7, t15, pkin(2) * t20 + 0; -t15 * t8, t15 * t7, -t18, pkin(2) * t22 + 0; -t7, -t8, 0, t19; 0, 0, 0, 1; t3, -t21, t15, t18 * pkin(1) + 0; t2, -t23, -t18, t15 * pkin(1) + 0; t13, t16, 0, t11; 0, 0, 0, 1; -t18 * t5, t18 * t4, t15, t18 * t1 + 0; -t15 * t5, t15 * t4, -t18, t15 * t1 + 0; -t4, -t5, 0, -pkin(3) * t7 + t19; 0, 0, 0, 1; t3, -t21, t15, t18 * t6 + 0; t2, -t23, -t18, t15 * t6 + 0; t13, t16, 0, t13 * pkin(4) + t11; 0, 0, 0, 1;];
T_ges = t10;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
