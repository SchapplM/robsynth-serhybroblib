% Calculate minimal parameter regressor of joint inertia matrix for
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
% MM_reg [((6+1)*6/2)x38]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 18:09
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = palh2m2OL_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 18:06:37
% EndTime: 2020-06-30 18:06:44
% DurationCPUTime: 0.84s
% Computational Cost: add. (816->89), mult. (1616->137), div. (0->0), fcn. (2130->10), ass. (0->68)
t44 = cos(qJ(2));
t31 = -t44 * pkin(4) - pkin(1);
t41 = sin(qJ(3));
t42 = sin(qJ(2));
t60 = cos(qJ(3));
t46 = t41 * t42 - t60 * t44;
t16 = t46 * pkin(2) + t31;
t23 = t41 * t44 + t60 * t42;
t40 = sin(qJ(4));
t59 = cos(qJ(4));
t45 = t40 * t23 + t59 * t46;
t7 = t45 * pkin(5) + t16;
t70 = 0.2e1 * t7;
t69 = 0.2e1 * t16;
t68 = 0.2e1 * t31;
t38 = sin(qJ(6));
t67 = -0.2e1 * t38;
t66 = 0.2e1 * t44;
t13 = t59 * t23 - t40 * t46;
t39 = sin(qJ(5));
t58 = cos(qJ(5));
t4 = t39 * t13 + t58 * t45;
t65 = t38 * t4;
t6 = t58 * t13 - t39 * t45;
t64 = t38 * t6;
t63 = t39 * pkin(5);
t62 = t41 * pkin(4);
t43 = cos(qJ(6));
t61 = t43 * t4;
t57 = t38 * t43;
t56 = pkin(3) * t38;
t33 = pkin(3) * t43;
t55 = t40 * pkin(2);
t54 = -0.2e1 * t6 * t4;
t53 = t6 * t57;
t35 = t60 * pkin(4);
t30 = t35 + pkin(2);
t51 = t59 * t62;
t22 = t40 * t30 + t51;
t52 = t58 * t22;
t21 = t59 * t30 - t40 * t62;
t18 = pkin(5) + t21;
t15 = t58 * t18;
t10 = -t39 * t22 + t15;
t50 = t58 * t55;
t11 = t39 * t18 + t52;
t9 = pkin(3) + t10;
t49 = t11 * t4 + t6 * t9;
t32 = t59 * pkin(2);
t29 = t32 + pkin(5);
t24 = t58 * t29;
t19 = -t39 * t55 + t24;
t17 = pkin(3) + t19;
t20 = t39 * t29 + t50;
t48 = t17 * t6 + t20 * t4;
t34 = t58 * pkin(5);
t28 = t34 + pkin(3);
t47 = t28 * t6 + t4 * t63;
t37 = t43 ^ 2;
t36 = t38 ^ 2;
t27 = 0.2e1 * t57;
t26 = t28 * t43;
t14 = t17 * t43;
t8 = t9 * t43;
t3 = t6 ^ 2;
t2 = t4 * pkin(3) + t7;
t1 = (t36 - t37) * t6;
t5 = [1, 0, 0, t42 ^ 2, t42 * t66, 0, 0, 0, pkin(1) * t66, -0.2e1 * pkin(1) * t42, t23 ^ 2, -0.2e1 * t23 * t46, 0, 0, 0, t46 * t68, t23 * t68, t13 ^ 2, -0.2e1 * t13 * t45, 0, 0, 0, t45 * t69, t13 * t69, t3, t54, 0, 0, 0, t4 * t70, t6 * t70, t37 * t3, -0.2e1 * t3 * t57, t43 * t54, 0.2e1 * t4 * t64, t4 ^ 2, 0.2e1 * t2 * t61, -0.2e1 * t2 * t65; 0, 0, 0, 0, 0, t42, t44, 0, 0, 0, 0, 0, t23, -t46, 0, 0, 0, 0, 0, t13, -t45, 0, 0, 0, 0, 0, t6, -t4, 0, 0, 0, -t53, t1, t65, t61, 0, t49 * t38, t49 * t43; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t62, 0, 0, 0, 0, 1, 0.2e1 * t21, -0.2e1 * t22, 0, 0, 0, 0, 1, 0.2e1 * t10, -0.2e1 * t11, t36, t27, 0, 0, 0, 0.2e1 * t8, t9 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t46, 0, 0, 0, 0, 0, t13, -t45, 0, 0, 0, 0, 0, t6, -t4, 0, 0, 0, -t53, t1, t65, t61, 0, t48 * t38, t48 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t62, 0, 0, 0, 0, 1, t21 + t32, -t51 + (-pkin(2) - t30) * t40, 0, 0, 0, 0, 1, t15 + t24 + (-t22 - t55) * t39, -t50 - t52 + (-t18 - t29) * t39, t36, t27, 0, 0, 0, t14 + t8, (-t17 - t9) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t32, -0.2e1 * t55, 0, 0, 0, 0, 1, 0.2e1 * t19, -0.2e1 * t20, t36, t27, 0, 0, 0, 0.2e1 * t14, t17 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t45, 0, 0, 0, 0, 0, t6, -t4, 0, 0, 0, -t53, t1, t65, t61, 0, t47 * t38, t47 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t21, -t22, 0, 0, 0, 0, 1, t10 + t34, -t52 + (-pkin(5) - t18) * t39, t36, t27, 0, 0, 0, t26 + t8, (-t28 - t9) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t32, -t55, 0, 0, 0, 0, 1, t19 + t34, -t50 + (-pkin(5) - t29) * t39, t36, t27, 0, 0, 0, t26 + t14, (-t17 - t28) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t34, -0.2e1 * t63, t36, t27, 0, 0, 0, 0.2e1 * t26, t28 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t4, 0, 0, 0, -t53, t1, t65, t61, 0, t6 * t56, t6 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t10, -t11, t36, t27, 0, 0, 0, t33 + t8, (-pkin(3) - t9) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t19, -t20, t36, t27, 0, 0, 0, t33 + t14, (-pkin(3) - t17) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t34, -t63, t36, t27, 0, 0, 0, t33 + t26, (-pkin(3) - t28) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t36, t27, 0, 0, 0, 0.2e1 * t33, -0.2e1 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 * t6, -t64, -t4, -t43 * t2, t38 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t43, 0, -t38 * t11, -t43 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t43, 0, -t38 * t20, -t43 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t43, 0, -t38 * t63, -t43 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t43, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
