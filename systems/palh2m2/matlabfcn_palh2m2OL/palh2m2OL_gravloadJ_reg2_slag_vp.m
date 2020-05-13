% Calculate inertial parameters regressor of gravitation load for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh2m2OL_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2OL_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t53 = sin(qJ(1));
t57 = cos(qJ(1));
t23 = g(1) * t57 + g(2) * t53;
t51 = sin(qJ(3));
t55 = cos(qJ(3));
t11 = g(3) * t51 + t23 * t55;
t52 = sin(qJ(2));
t56 = cos(qJ(2));
t92 = g(3) * t55 - t51 * t23;
t2 = t11 * t52 - t56 * t92;
t47 = qJ(2) + qJ(3);
t42 = qJ(4) + t47;
t35 = sin(t42);
t69 = t23 * t35;
t36 = cos(t42);
t76 = g(3) * t36;
t6 = t69 - t76;
t91 = pkin(5) * t6;
t37 = qJ(5) + t42;
t30 = sin(t37);
t90 = t23 * t30;
t31 = cos(t37);
t87 = g(3) * t30 + t23 * t31;
t41 = cos(t47);
t86 = pkin(2) * (t23 * sin(t47) - t41 * g(3));
t9 = -g(3) * t56 + t23 * t52;
t84 = pkin(4) * t9;
t29 = pkin(5) * t36;
t77 = g(3) * t31;
t4 = -t77 + t90;
t50 = sin(qJ(6));
t72 = t50 * t4;
t45 = t56 * pkin(4);
t71 = pkin(5) * t69 + t86;
t24 = -g(1) * t53 + g(2) * t57;
t67 = t50 * t24;
t54 = cos(qJ(6));
t66 = t54 * t24;
t34 = pkin(2) * t41;
t63 = t31 * pkin(3);
t62 = t45 + pkin(1);
t61 = t34 + t62;
t60 = t29 + t61;
t38 = qJ(1) + t42;
t39 = qJ(1) - t42;
t43 = qJ(1) + t47;
t44 = qJ(1) - t47;
t48 = qJ(2) + qJ(1);
t49 = -qJ(2) + qJ(1);
t58 = (-t34 - t45 - t29) * g(3) + ((-cos(t48) / 0.2e1 + cos(t49) / 0.2e1) * pkin(4) + (-cos(t38) / 0.2e1 + cos(t39) / 0.2e1) * pkin(5) + (-cos(t43) / 0.2e1 + cos(t44) / 0.2e1) * pkin(2)) * g(2) + ((sin(t38) / 0.2e1 - sin(t39) / 0.2e1) * pkin(5) + (sin(t48) / 0.2e1 - sin(t49) / 0.2e1) * pkin(4) + (sin(t43) / 0.2e1 - sin(t44) / 0.2e1) * pkin(2)) * g(1);
t33 = qJ(1) - t37;
t32 = qJ(1) + t37;
t13 = t30 * t24;
t10 = pkin(3) * t90;
t7 = g(3) * t35 + t23 * t36;
t3 = t54 * t4;
t1 = t11 * t56 + t52 * t92;
t5 = [0, 0, 0, 0, 0, 0, -t24, t23, 0, 0, 0, 0, 0, 0, 0, 0, -t56 * t24, t52 * t24, -t23, -pkin(1) * t24, 0, 0, 0, 0, 0, 0, t24 * (t51 * t52 - t55 * t56), t24 * (t51 * t56 + t55 * t52), -t23, -t24 * t62, 0, 0, 0, 0, 0, 0, -t36 * t24, t35 * t24, -t23, -t61 * t24, 0, 0, 0, 0, 0, 0, -t31 * t24, t13, -t23, -t60 * t24, 0, 0, 0, 0, 0, 0, t50 * t23 - t31 * t66, t54 * t23 + t31 * t67, t13, -t24 * (t60 + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, g(3) * t52 + t23 * t56, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t84, 0, 0, 0, 0, 0, 0, t6, t7, 0, t2 * pkin(2) + t84, 0, 0, 0, 0, 0, 0, t4, t87, 0, t58, 0, 0, 0, 0, 0, 0, t3, -t72, t87, t58 + (-t77 + (cos(t33) / 0.2e1 - cos(t32) / 0.2e1) * g(2) + (sin(t32) / 0.2e1 - sin(t33) / 0.2e1) * g(1)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, t86, 0, 0, 0, 0, 0, 0, t4, t87, 0, -pkin(5) * t76 + t71, 0, 0, 0, 0, 0, 0, t3, -t72, t87, t10 + (-t63 - t29) * g(3) + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, 0, 0, 0, t4, t87, 0, t91, 0, 0, 0, 0, 0, 0, t3, -t72, t87, -g(3) * t63 + t10 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t87, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t72, t87, pkin(3) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t87 - t66, t54 * t87 + t67, 0, 0;];
taug_reg = t5;
