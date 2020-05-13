% Calculate inertial parameters regressor of gravitation load for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% taug_reg [10x(10*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh3m2OL_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_gravloadJ_reg2_slag_vp: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_gravloadJ_reg2_slag_vp: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t46 = sin(qJ(1));
t50 = cos(qJ(1));
t74 = -g(1) * t46 + g(2) * t50;
t24 = g(1) * t50 + g(2) * t46;
t42 = qJ(2) + qJ(3);
t39 = qJ(4) + t42;
t31 = sin(t39);
t32 = cos(t39);
t73 = -g(3) * t32 + t24 * t31;
t35 = sin(t42);
t72 = pkin(4) * t35;
t37 = cos(t42);
t71 = pkin(4) * t37;
t70 = pkin(8) * t31;
t69 = pkin(10) * t32;
t66 = g(3) * t31;
t45 = sin(qJ(2));
t64 = t45 * pkin(1);
t44 = sin(qJ(5));
t63 = t46 * t44;
t48 = cos(qJ(5));
t62 = t46 * t48;
t61 = t50 * t44;
t60 = t50 * t48;
t49 = cos(qJ(2));
t40 = t49 * pkin(1);
t59 = t40 + pkin(12);
t41 = qJ(2) + qJ(7);
t58 = t40 - t71;
t21 = -t64 + t72;
t57 = t21 - t69;
t38 = pkin(15) - t41;
t56 = -t69 + t72;
t55 = t32 * pkin(8) + t31 * pkin(10);
t52 = -t55 - t71;
t8 = g(3) * t37 - t24 * t35;
t51 = -g(3) * t49 + t24 * t45;
t47 = cos(qJ(6));
t43 = sin(qJ(6));
t36 = cos(t41);
t34 = sin(t41);
t33 = -qJ(8) + t38;
t29 = cos(t38);
t28 = sin(t38);
t27 = cos(t33);
t26 = sin(t33);
t25 = pkin(3) * t29;
t23 = t50 * t70;
t22 = t46 * t70;
t20 = pkin(12) + t58;
t17 = t32 * t60 - t63;
t16 = t32 * t61 + t62;
t15 = t32 * t62 + t61;
t14 = t32 * t63 - t60;
t13 = t74 * t31;
t12 = t74 * t59;
t11 = t51 * pkin(1);
t10 = -g(3) * t35 - t24 * t37;
t9 = g(3) * t34 + t24 * t36;
t7 = -g(3) * t36 + t24 * t34;
t6 = t24 * t32 + t66;
t4 = g(3) * t26 - t24 * t27;
t3 = g(3) * t27 + t24 * t26;
t2 = t73 * t48;
t1 = t73 * t44;
t5 = [0, 0, 0, 0, 0, 0, -t74, t24, 0, 0, 0, 0, 0, 0, 0, 0, -t74 * t49, t74 * t45, -t24, -t74 * pkin(12), 0, 0, 0, 0, 0, 0, t74 * t37, -t74 * t35, -t24, -t12, 0, 0, 0, 0, 0, 0, t74 * t32, -t13, -t24, -t74 * t20, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t17, g(1) * t14 - g(2) * t16, t13, t74 * (-t20 + t55), 0, 0, 0, 0, 0, 0, -t74 * t47, t74 * t43, -t24, t74 * pkin(6), 0, 0, 0, 0, 0, 0, -t74 * t36, t74 * t34, -t24, -t12, 0, 0, 0, 0, 0, 0, t74 * t27, t74 * t26, -t24, -t74 * (t25 + t59); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, g(3) * t45 + t24 * t49, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10, 0, t11, 0, 0, 0, 0, 0, 0, -t73, -t6, 0, -g(3) * t58 - t24 * t21, 0, 0, 0, 0, 0, 0, -t2, t1, t6, -g(1) * (t57 * t50 + t23) - g(2) * (t57 * t46 + t22) - g(3) * (t40 + t52), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t9, 0, t11, 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * (t25 + t40) - t24 * (pkin(3) * t28 - t64); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t6, 0, t8 * pkin(4), 0, 0, 0, 0, 0, 0, -t2, t1, t6, -g(1) * (t56 * t50 + t23) - g(2) * (t56 * t46 + t22) - g(3) * t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t6, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, t6, -g(1) * (-t50 * t69 + t23) - g(2) * (-t46 * t69 + t22) + g(3) * t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14 - t44 * t66, -g(1) * t17 - g(2) * t15 - t48 * t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t47 + t24 * t43, g(3) * t43 + t24 * t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t9, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, (-g(3) * t29 - t24 * t28) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t5;
