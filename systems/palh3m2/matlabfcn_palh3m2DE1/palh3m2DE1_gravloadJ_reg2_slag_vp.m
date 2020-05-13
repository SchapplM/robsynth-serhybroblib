% Calculate inertial parameters regressor of gravitation load for
% palh3m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh3m2DE1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_gravloadJ_reg2_slag_vp: pkin has to be [18x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t69 = sin(qJ(1));
t64 = g(2) * t69;
t55 = cos(qJ(1));
t73 = g(1) * t55;
t27 = t64 + t73;
t48 = sin(qJ(3));
t53 = cos(qJ(3));
t14 = t53 * g(3) - t48 * t27;
t54 = cos(qJ(2));
t76 = t14 * t54;
t71 = t48 * g(3);
t15 = t27 * t53 + t71;
t49 = sin(qJ(2));
t4 = -t49 * t15 + t76;
t74 = t4 * pkin(4);
t72 = g(3) * t54;
t68 = t48 * t49;
t28 = -g(1) * t69 + g(2) * t55;
t52 = cos(qJ(4));
t67 = t52 * t28;
t66 = pkin(1) * t54 + pkin(12);
t32 = -pkin(4) * t53 + pkin(1);
t65 = pkin(4) * t68 + t32 * t54 + pkin(12);
t16 = g(3) * t49 + t27 * t54;
t50 = sin(pkin(15));
t56 = cos(pkin(15));
t17 = -t50 * g(3) + t27 * t56;
t19 = g(3) * t56 + t50 * t27;
t40 = pkin(17) + pkin(18);
t33 = sin(t40);
t34 = cos(t40);
t42 = sin(pkin(16));
t43 = cos(pkin(16));
t63 = -(t17 * t43 - t19 * t42) * t34 + (t17 * t42 + t19 * t43) * t33;
t62 = t49 * t27 - t72;
t21 = t42 * t56 + t43 * t50;
t22 = -t42 * t50 + t43 * t56;
t61 = t28 * (t34 * t21 + t33 * t22);
t60 = t33 * t21 - t22 * t34;
t44 = sin(pkin(18));
t46 = cos(pkin(18));
t59 = t44 * t56 + t46 * t50;
t58 = t60 * t28;
t57 = cos(pkin(14));
t51 = sin(pkin(14));
t47 = sin(qJ(4));
t41 = qJ(3) + qJ(2);
t36 = cos(t41);
t35 = sin(t41);
t26 = t56 * pkin(8) - t50 * pkin(10);
t25 = t50 * pkin(8) + t56 * pkin(10);
t23 = -t44 * t50 + t46 * t56;
t20 = t66 * t28;
t11 = pkin(1) * t62;
t10 = -g(3) * t35 - t27 * t36;
t9 = g(3) * t36 - t27 * t35;
t8 = t16 * t50 + t56 * t62;
t7 = -t49 * t14 - t15 * t54;
t6 = t16 * t56 - t50 * t62;
t1 = (pkin(1) * t64 + t32 * t73) * t49 - pkin(1) * t72 + ((-t53 * t64 - t71) * t49 + t76) * pkin(4);
t2 = [0, 0, 0, 0, 0, 0, -t28, t27, 0, 0, 0, 0, 0, 0, 0, 0, -t54 * t28, t49 * t28, -t27, -pkin(12) * t28, 0, 0, 0, 0, 0, 0, -t28 * (-t53 * t54 + t68), -t28 * (t48 * t54 + t53 * t49), -t27, -t20, 0, 0, 0, 0, 0, 0, -t58, t61, -t27, -t28 * t65, 0, 0, 0, 0, 0, 0, -t47 * t27 - t60 * t67, -t52 * t27 + t47 * t58, -t61, -t28 * ((t42 * t25 - t26 * t43) * t34 + (t25 * t43 + t42 * t26) * t33 + t65), 0, 0, 0, 0, 0, 0, -t28 * ((-t50 * t49 + t54 * t56) * t57 + t51 * (t49 * t56 + t54 * t50)), -t28 * ((-t51 * t50 - t57 * t56) * t49 + (-t57 * t50 + t51 * t56) * t54), -t27, pkin(6) * t28, 0, 0, 0, 0, 0, 0, t23 * t28, t59 * t28, -t27, -t20, 0, 0, 0, 0, 0, 0, t36 * t28, -t35 * t28, -t27, -t28 * ((-cos(pkin(17)) * t23 + t59 * sin(pkin(17))) * pkin(3) + t66); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t16, 0, 0, 0, 0, 0, 0, 0, 0, t4, t7, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, -t51 * t6 + t8 * t57, t51 * t8 + t6 * t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, t9, t10, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t47 + t67, -t47 * t28 + t63 * t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t2;
