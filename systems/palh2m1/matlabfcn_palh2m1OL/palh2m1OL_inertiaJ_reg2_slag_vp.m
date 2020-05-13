% Calculate inertial parameters regressor of joint inertia matrix for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = palh2m1OL_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t48 = cos(qJ(3));
t38 = t48 * pkin(2);
t33 = t38 + pkin(3);
t43 = sin(qJ(4));
t47 = cos(qJ(4));
t44 = sin(qJ(3));
t70 = t44 * pkin(2);
t58 = t47 * t70;
t18 = t33 * t43 + t58;
t16 = pkin(6) + t18;
t42 = sin(qJ(5));
t40 = t42 ^ 2;
t46 = cos(qJ(5));
t41 = t46 ^ 2;
t61 = t40 + t41;
t64 = t61 * t16;
t71 = t43 * pkin(3);
t31 = pkin(6) + t71;
t78 = t61 * t31;
t49 = cos(qJ(2));
t45 = sin(qJ(2));
t66 = t44 * t45;
t21 = -t48 * t49 + t66;
t22 = t44 * t49 + t48 * t45;
t10 = t43 * t21 - t47 * t22;
t77 = -0.2e1 * t10;
t8 = -t47 * t21 - t43 * t22;
t76 = t8 ^ 2;
t15 = -pkin(3) * t66 + pkin(1) + (t48 * pkin(3) + pkin(2)) * t49;
t75 = 0.2e1 * t15;
t34 = t49 * pkin(2) + pkin(1);
t74 = -0.2e1 * t34;
t73 = 0.2e1 * t49;
t72 = pkin(4) * t42;
t5 = t42 * t8;
t6 = t46 * t8;
t56 = -t33 * t47 + t43 * t70;
t17 = -pkin(4) + t56;
t69 = t17 * t46;
t37 = t47 * pkin(3);
t32 = -t37 - pkin(4);
t68 = t32 * t46;
t67 = t42 * t46;
t65 = t46 * t10;
t62 = pkin(6) * t61;
t60 = -0.2e1 * t5;
t59 = 0.2e1 * t6;
t55 = -pkin(4) * t10 - pkin(6) * t8;
t54 = t10 * t17 - t16 * t8;
t53 = t10 * t32 - t31 * t8;
t39 = pkin(4) * t46;
t29 = 0.2e1 * t67;
t26 = t32 * t42;
t13 = t17 * t42;
t7 = t10 ^ 2;
t4 = t42 * t65;
t3 = (-t40 + t41) * t10;
t2 = t8 * pkin(4) - t10 * pkin(6) + t15;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t45 ^ 2, t45 * t73, 0, t49 ^ 2, 0, 0, pkin(1) * t73, -0.2e1 * pkin(1) * t45, 0, pkin(1) ^ 2, t22 ^ 2, -0.2e1 * t22 * t21, 0, t21 ^ 2, 0, 0, t21 * t74, t22 * t74, 0, t34 ^ 2, t7, t8 * t77, 0, t76, 0, 0, t8 * t75, t10 * t75, 0, t15 ^ 2, t41 * t7, -0.2e1 * t7 * t67, t10 * t59, t40 * t7, t10 * t60, t76, t2 * t59, t2 * t60, t61 * t2 * t77, t61 * t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, -t49, 0, 0, 0, 0, 0, 0, 0, -t22, 0, t21, 0, 0, 0, (t21 * t44 + t22 * t48) * pkin(2), 0, 0, 0, t10, 0, -t8, 0, 0, 0, t10 * t56 - t18 * t8, 0, t4, t3, t5, -t4, t6, 0, t54 * t42, t54 * t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t70, 0, (t44 ^ 2 + t48 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t56, -0.2e1 * t18, 0, t18 ^ 2 + t56 ^ 2, t40, t29, 0, t41, 0, 0, -0.2e1 * t69, 0.2e1 * t13, 0.2e1 * t64, t61 * t16 ^ 2 + t17 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, t21, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t8, 0, 0, 0, (-t10 * t47 - t43 * t8) * pkin(3), 0, t4, t3, t5, -t4, t6, 0, t53 * t42, t53 * t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t38, -t70, 0, 0, 0, 0, 0, 0, 0, 1, t37 - t56, -t58 + (-pkin(3) - t33) * t43, 0, (t18 * t43 - t47 * t56) * pkin(3), t40, t29, 0, t41, 0, 0, (-t17 - t32) * t46, t26 + t13, t78 + t64, t16 * t78 + t17 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t37, -0.2e1 * t71, 0, (t43 ^ 2 + t47 ^ 2) * pkin(3) ^ 2, t40, t29, 0, t41, 0, 0, -0.2e1 * t68, 0.2e1 * t26, 0.2e1 * t78, t61 * t31 ^ 2 + t32 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t8, 0, 0, 0, 0, 0, t4, t3, t5, -t4, t6, 0, t55 * t42, t55 * t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t56, -t18, 0, 0, t40, t29, 0, t41, 0, 0, t39 - t69, t13 - t72, t62 + t64, -t17 * pkin(4) + pkin(6) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t37, -t71, 0, 0, t40, t29, 0, t41, 0, 0, t39 - t68, t26 - t72, t62 + t78, -t32 * pkin(4) + pkin(6) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t40, t29, 0, t41, 0, 0, 0.2e1 * t39, -0.2e1 * t72, 0.2e1 * t62, t61 * pkin(6) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, -t42 * t10, t8, t46 * t2, -t42 * t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t46, 0, -t42 * t16, -t46 * t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t46, 0, -t42 * t31, -t46 * t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t46, 0, -t42 * pkin(6), -t46 * pkin(6), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
