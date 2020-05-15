% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% fourbar1DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
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
% Datum: 2020-04-24 20:05
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = fourbar1DE2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:02:24
% EndTime: 2020-04-24 20:02:25
% DurationCPUTime: 0.25s
% Computational Cost: add. (644->48), mult. (847->54), div. (66->3), fcn. (253->10), ass. (0->48)
t27 = 0.1e1 / pkin(4);
t24 = cos(qJ(1));
t50 = pkin(2) * t24;
t40 = pkin(1) * t50;
t21 = -0.2e1 * t40;
t30 = pkin(2) ^ 2;
t31 = pkin(1) ^ 2;
t41 = t30 + t31;
t18 = t21 + t41;
t17 = 0.1e1 / t18;
t47 = t17 / 0.2e1;
t35 = t27 * t47;
t26 = pkin(4) ^ 2;
t28 = pkin(3) ^ 2;
t36 = -t28 + t41;
t16 = t21 + t26 + t36;
t23 = sin(qJ(1));
t22 = t23 * pkin(2);
t42 = t21 + t31;
t52 = -pkin(3) + pkin(4);
t53 = -pkin(3) - pkin(4);
t13 = sqrt(-((pkin(2) - t53) * (pkin(2) + t53) + t42) * ((pkin(2) - t52) * (pkin(2) + t52) + t42));
t19 = -pkin(1) + t50;
t45 = t19 * t13;
t55 = -t16 * t22 + t45;
t57 = t55 * t35;
t29 = 0.1e1 / pkin(3);
t46 = t17 * t29;
t34 = t46 / 0.2e1;
t37 = -t26 + t28 + t31;
t33 = t30 + t37;
t15 = t21 + t33;
t48 = t15 * t23;
t56 = (pkin(2) * t48 + t45) * t34;
t51 = pkin(1) * t24;
t49 = t13 * t23;
t44 = t27 * t29;
t20 = -pkin(2) + t51;
t4 = qJ(1) + atan2((pkin(1) * t48 - t20 * t13) * t46, (-pkin(1) * t49 - t20 * t15) * t46);
t38 = t22 + 0;
t32 = 0 + t50;
t12 = pkin(2) * t49;
t7 = (t19 * t16 + t12) * t35;
t6 = (-t19 * t15 + t12) * t34;
t3 = atan2(t13 * t44, (t26 - t36 + 0.2e1 * t40) * t44) + t4;
t2 = cos(t3);
t1 = sin(t3);
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t24, -t23, 0, 0; t23, t24, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t6, t56, 0, t32; -t56, t6, 0, t38; 0, 0, 1, 0; 0, 0, 0, 1; t7, t57, 0, pkin(1) + 0; -t57, t7, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t2, t1, 0, -pkin(3) * cos(t4) + t32; -t1, -t2, 0, -pkin(3) * sin(t4) + t38; 0, 0, 1, 0; 0, 0, 0, 1; t7, t57, 0, (t12 + 0.2e1 * t41 * 0 + pkin(1) * t33 + (-pkin(2) * (0.4e1 * pkin(1) * 0 - t30 + t37) - 0.2e1 * t30 * t51) * t24) * t47; -t57, t7, 0, (0.2e1 * 0 * t18 - t55) * t47; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t5;
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
