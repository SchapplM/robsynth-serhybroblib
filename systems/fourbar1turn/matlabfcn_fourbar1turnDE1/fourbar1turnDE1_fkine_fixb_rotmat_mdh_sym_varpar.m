% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% fourbar1turnDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
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
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = fourbar1turnDE1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:25:20
% EndTime: 2020-04-12 19:25:20
% DurationCPUTime: 0.74s
% Computational Cost: add. (4021->56), mult. (5771->73), div. (411->3), fcn. (1690->15), ass. (0->57)
t62 = -pkin(3) - pkin(4);
t61 = -pkin(3) + pkin(4);
t33 = cos(qJ(2));
t60 = pkin(2) * t33;
t31 = sin(qJ(2));
t27 = t31 * pkin(2);
t32 = sin(qJ(1));
t49 = pkin(1) * t60;
t26 = -0.2e1 * t49;
t40 = pkin(1) ^ 2;
t51 = t26 + t40;
t18 = sqrt(-((pkin(2) - t62) * (pkin(2) + t62) + t51) * ((pkin(2) - t61) * (pkin(2) + t61) + t51));
t35 = pkin(4) ^ 2;
t37 = pkin(3) ^ 2;
t50 = pkin(2) ^ 2 + t40;
t43 = -t37 + t50;
t20 = t26 + t35 + t43;
t24 = pkin(1) - t60;
t44 = t26 + t50;
t21 = 0.1e1 / t44;
t36 = 0.1e1 / pkin(4);
t56 = t21 * t36;
t57 = t18 * t31;
t11 = atan2((t24 * t18 + t20 * t27) * t56 / 0.2e1, -(-pkin(2) * t57 + t24 * t20) * t56 / 0.2e1);
t9 = sin(t11);
t59 = t32 * t9;
t34 = cos(qJ(1));
t58 = t34 * t9;
t38 = 0.1e1 / pkin(3);
t55 = t21 * t38;
t10 = cos(t11);
t7 = t32 * t10;
t54 = t32 * t33;
t8 = t34 * t10;
t53 = t34 * t33;
t52 = t36 * t38;
t30 = pkin(5) + 0;
t48 = pkin(2) * t54 + 0;
t47 = pkin(2) * t53 + 0;
t46 = t32 * pkin(1) + 0;
t45 = t34 * pkin(1) + 0;
t42 = t27 + t30;
t19 = -t35 + t37 + t44;
t25 = pkin(1) * t33 - pkin(2);
t14 = atan2((pkin(1) * t31 * t19 - t25 * t18) * t55, (-pkin(1) * t57 - t25 * t19) * t55);
t12 = sin(t14);
t13 = cos(t14);
t41 = t33 * t12 + t31 * t13;
t5 = t31 * t12 - t33 * t13;
t17 = atan2(t18 * t52, (t35 - t43 + 0.2e1 * t49) * t52);
t16 = cos(t17);
t15 = sin(t17);
t4 = t5 * t34;
t3 = t41 * t34;
t2 = t5 * t32;
t1 = t41 * t32;
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t32, 0, 0; t32, t34, 0, 0; 0, 0, 1, t30; 0, 0, 0, 1; t53, -t34 * t31, t32, 0; t54, -t32 * t31, -t34, 0; t31, t33, 0, t30; 0, 0, 0, 1; t4, t3, t32, t47; t2, t1, -t34, t48; -t41, t5, 0, t42; 0, 0, 0, 1; t8, -t58, t32, t45; t7, -t59, -t34, t46; t9, t10, 0, t30; 0, 0, 0, 1; t3 * t15 + t4 * t16, -t4 * t15 + t3 * t16, t32, t4 * pkin(3) + t47; t1 * t15 + t2 * t16, t1 * t16 - t2 * t15, -t34, t2 * pkin(3) + t48; t5 * t15 - t16 * t41, t15 * t41 + t5 * t16, 0, -pkin(3) * t41 + t42; 0, 0, 0, 1; t8, -t58, t32, pkin(4) * t8 + t45; t7, -t59, -t34, pkin(4) * t7 + t46; t9, t10, 0, t9 * pkin(4) + t30; 0, 0, 0, 1;];
T_ges = t6;
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
