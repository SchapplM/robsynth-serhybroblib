% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% fourbarprisDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   4:  mdh base (link 0) -> mdh frame (4-1), link (4-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = fourbarprisDE2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:43:52
% EndTime: 2020-05-07 09:43:53
% DurationCPUTime: 0.18s
% Computational Cost: add. (360->26), mult. (240->42), div. (72->3), fcn. (50->6), ass. (0->39)
t19 = (qJ(1) ^ 2);
t42 = (-pkin(3) ^ 2 - t19);
t23 = (pkin(1) ^ 2);
t41 = (t23 + t42);
t40 = 2 * pkin(1);
t24 = 1 / pkin(1);
t39 = t24 / 0.2e1;
t38 = (pkin(3) * qJ(1));
t18 = (qJ(1) + pkin(3));
t17 = 1 / t18;
t22 = 1 / pkin(2);
t37 = t17 * t22;
t36 = t17 * t24;
t35 = t22 * t24;
t33 = (-pkin(2) - t18);
t32 = (-pkin(2) + t18);
t31 = -t36 / 0.2e1;
t30 = t36 / 0.2e1;
t29 = -t35 / 0.2e1;
t28 = t35 / 0.2e1;
t27 = 0 * t40;
t26 = -2 * t38 + t41;
t21 = (pkin(2) ^ 2);
t12 = t21 + t26;
t15 = pkin(1) + 0;
t14 = 0 * t40;
t13 = -t21 + t23 + 2 * t38 - t42;
t11 = t12 * t28;
t10 = t13 * t31;
t9 = (t14 + t12) * t39;
t8 = sqrt(-((pkin(1) + t33) * (pkin(1) + t32) * (pkin(1) - t32) * (pkin(1) - t33)));
t7 = t8 * t28;
t6 = t8 * t31;
t5 = t8 * t30;
t4 = (-t8 + t27) * t39;
t3 = atan2(t7, t12 * t29) + atan2(t8 * t37, (-t21 + t26) * t37);
t2 = cos(t3);
t1 = sin(t3);
t16 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t10, t5, 0, t15; t6, t10, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t6, 0, t10, (2 * (pkin(1) * t15 - t19) * pkin(3) + (t14 + t21 + t41) * qJ(1)) * t30; t13 * t30, 0, t6, (-t8 * qJ(1) + (t18 * t27)) * t30; 0, -1, 0, 0; 0, 0, 0, 1; t11, t7, 0, 0; t8 * t29, t11, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t2, -t1, 0, t9; t1, t2, 0, t4; 0, 0, 1, 0; 0, 0, 0, 1; t10, t5, 0, t9; t6, t10, 0, t4; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t16;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
