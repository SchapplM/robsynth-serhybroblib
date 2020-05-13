% Calculate inertial parameters regressor of gravitation load for
% palh1m2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh1m2DE1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_gravloadJ_reg2_slag_vp: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t67 = cos(qJ(1));
t82 = sin(qJ(1));
t45 = g(1) * t67 + g(2) * t82;
t62 = sin(pkin(18));
t68 = cos(pkin(18));
t26 = t62 * g(3) + t45 * t68;
t28 = g(3) * t68 - t62 * t45;
t51 = sin(pkin(22));
t55 = cos(pkin(22));
t85 = -t26 * t55 + t51 * t28;
t32 = t62 * t51 + t55 * t68;
t33 = -t68 * t51 + t62 * t55;
t52 = sin(pkin(21));
t56 = cos(pkin(21));
t15 = t32 * t56 + t52 * t33;
t16 = -t52 * t32 + t33 * t56;
t46 = -g(1) * t82 + g(2) * t67;
t53 = sin(pkin(20));
t57 = cos(pkin(20));
t84 = t46 * (-t53 * t15 + t16 * t57);
t65 = cos(qJ(3));
t83 = t65 * g(3);
t80 = t45 * t65;
t60 = sin(qJ(3));
t61 = sin(qJ(2));
t66 = cos(qJ(2));
t75 = t65 * t66;
t70 = t61 * t60 - t75;
t79 = t46 * t70;
t59 = sin(qJ(4));
t77 = t59 * t46;
t64 = cos(qJ(4));
t76 = t64 * t46;
t74 = pkin(1) * t61 - pkin(15);
t23 = -t60 * t45 + t83;
t24 = t60 * g(3) + t80;
t11 = t61 * t23 + t66 * t24;
t12 = -t23 * t66 + t24 * t61;
t54 = sin(pkin(19));
t58 = cos(pkin(19));
t4 = t54 * t11 + t12 * t58;
t73 = g(3) * t61 + t45 * t66;
t72 = g(3) * t66 - t61 * t45;
t43 = -t62 * pkin(9) + t68 * pkin(11);
t44 = t68 * pkin(9) + t62 * pkin(11);
t71 = t51 * t43 - t44 * t55;
t69 = cos(pkin(17));
t63 = sin(pkin(17));
t47 = pkin(5) * t60 + pkin(1);
t42 = pkin(15) * t46;
t41 = t66 * t46;
t39 = t61 * t46;
t37 = t63 * t62 + t69 * t68;
t36 = t66 * t60 + t65 * t61;
t35 = t62 * t69 - t68 * t63;
t31 = t54 * t66 + t61 * t58;
t30 = t61 * t54 - t58 * t66;
t29 = t74 * t46;
t20 = pkin(1) * t73;
t18 = t46 * t36;
t17 = t43 * t55 + t51 * t44;
t14 = t62 * t73 + t68 * t72;
t13 = -t62 * t72 + t68 * t73;
t10 = t51 * t26 + t28 * t55;
t9 = t12 * pkin(5);
t7 = (-pkin(5) * t83 + t45 * t47) * t66 + t61 * (pkin(5) * t80 + g(3) * t47);
t5 = t15 * t57 + t53 * t16;
t3 = t11 * t58 - t54 * t12;
t2 = t4 * pkin(2);
t1 = (t10 * t52 + t85 * t56) * t57 + (t10 * t56 - t85 * t52) * t53;
t6 = [0, 0, 0, 0, 0, 0, -t46, t45, 0, 0, 0, 0, 0, 0, 0, 0, t39, t41, -t45, -t42, 0, 0, 0, 0, 0, 0, t79, t18, -t45, t29, 0, 0, 0, 0, 0, 0, t46 * t5, -t84, -t45, t46 * (t70 * pkin(5) + t74), 0, 0, 0, 0, 0, 0, -t59 * t45 + t5 * t76, -t64 * t45 - t5 * t77, t84, ((-t17 * t52 - t71 * t56) * t57 + t61 * t47 - pkin(5) * t75 - pkin(15) + (-t17 * t56 + t71 * t52) * t53) * t46, 0, 0, 0, 0, 0, 0, (-t35 * t66 + t61 * t37) * t46, t46 * (t61 * t35 + t37 * t66), -t45, pkin(14) * t46, 0, 0, 0, 0, 0, 0, t32 * t46, -t33 * t46, -t45, t29, 0, 0, 0, 0, 0, 0, t46 * (t30 * t65 + t60 * t31), -t46 * (t60 * t30 - t65 * t31), -t45, -t42, 0, 0, 0, 0, 0, 0, t39, t41, -t45, -(pkin(15) + (-t36 * t54 - t58 * t70) * pkin(2)) * t46, 0, 0, 0, 0, 0, 0, t79, t18, -t45, t46 * (t15 * pkin(4) + t74); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t72, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, t13 * t69 + t14 * t63, -t13 * t63 + t14 * t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, 0, 0, 0, t73, t72, 0, t2, 0, 0, 0, 0, 0, 0, t12, t11, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, t12, t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t59 + t76, t1 * t64 - t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t6;
