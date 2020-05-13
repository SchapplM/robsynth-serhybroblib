% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% fourbar1turnTE
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
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = fourbar1turnTE_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:18:25
% EndTime: 2020-04-12 19:18:25
% DurationCPUTime: 0.26s
% Computational Cost: add. (1356->52), mult. (1892->84), div. (135->3), fcn. (604->6), ass. (0->64)
t35 = pkin(4) ^ 2;
t37 = pkin(3) ^ 2;
t33 = cos(qJ(2));
t64 = pkin(2) * t33;
t56 = pkin(1) * t64;
t26 = -0.2e1 * t56;
t40 = pkin(1) ^ 2;
t57 = pkin(2) ^ 2 + t40;
t51 = t26 + t57;
t19 = -t35 + t37 + t51;
t25 = pkin(1) * t33 - pkin(2);
t58 = t26 + t40;
t65 = -pkin(3) + pkin(4);
t66 = -pkin(3) - pkin(4);
t17 = sqrt(-((pkin(2) - t66) * (pkin(2) + t66) + t58) * ((pkin(2) - t65) * (pkin(2) + t65) + t58));
t31 = sin(qJ(2));
t63 = t17 * t31;
t13 = -pkin(1) * t63 - t25 * t19;
t16 = pkin(1) * t31 * t19 - t25 * t17;
t21 = 0.1e1 / t51;
t38 = 0.1e1 / pkin(3);
t62 = t21 * t38;
t67 = t31 / 0.2e1;
t72 = (t33 * t16 / 0.2e1 + t13 * t67) * t62;
t71 = -t17 / 0.2e1;
t70 = t17 / 0.2e1;
t50 = -t37 + t57;
t69 = t35 / 0.2e1 - t50 / 0.2e1 + t56;
t68 = -t21 / 0.2e1;
t27 = t31 * pkin(2);
t32 = sin(qJ(1));
t61 = t32 * t33;
t34 = cos(qJ(1));
t60 = t34 * t33;
t36 = 0.1e1 / pkin(4);
t59 = t36 * t38;
t30 = pkin(5) + 0;
t55 = pkin(2) * t61 + 0;
t54 = pkin(2) * t60 + 0;
t53 = t32 * pkin(1) + 0;
t52 = t34 * pkin(1) + 0;
t20 = t26 + t35 + t50;
t24 = pkin(1) - t64;
t14 = -pkin(2) * t63 + t24 * t20;
t49 = t14 * t68;
t15 = t24 * t17 + t20 * t27;
t48 = t15 * t21 / 0.2e1;
t47 = t36 * t68;
t46 = t27 + t30;
t45 = t15 * t47;
t44 = t32 * t49;
t43 = t34 * t49;
t6 = (t16 * t67 - t33 * t13 / 0.2e1) * t62;
t12 = t36 * t48;
t11 = t14 * t47;
t10 = t34 * t45;
t9 = t36 * t43;
t8 = t32 * t45;
t7 = t36 * t44;
t4 = t34 * t72;
t3 = t34 * t6;
t2 = t32 * t72;
t1 = t32 * t6;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t32, 0, 0; t32, t34, 0, 0; 0, 0, 1, t30; 0, 0, 0, 1; t60, -t34 * t31, t32, 0; t61, -t32 * t31, -t34, 0; t31, t33, 0, t30; 0, 0, 0, 1; t3, t4, t32, t54; t1, t2, -t34, t55; -t72, t6, 0, t46; 0, 0, 0, 1; t9, t10, t32, t52; t7, t8, -t34, t53; t12, t11, 0, t30; 0, 0, 0, 1; (t3 * t69 + t4 * t70) * t59, (t3 * t71 + t4 * t69) * t59, t32, t3 * pkin(3) + t54; (t1 * t69 + t2 * t70) * t59, (t1 * t71 + t2 * t69) * t59, -t34, t1 * pkin(3) + t55; (t6 * t70 - t69 * t72) * t59, (t6 * t69 - t71 * t72) * t59, 0, -pkin(3) * t72 + t46; 0, 0, 0, 1; t9, t10, t32, t43 + t52; t7, t8, -t34, t44 + t53; t12, t11, 0, t48 + t30; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
