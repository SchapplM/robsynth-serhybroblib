% Calculate inertial parameters regressor of gravitation load for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh2m1OL_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t41 = sin(qJ(2));
t45 = cos(qJ(2));
t56 = sin(qJ(1));
t51 = g(2) * t56;
t46 = cos(qJ(1));
t62 = g(1) * t46;
t32 = t51 + t62;
t39 = sin(qJ(4));
t43 = cos(qJ(4));
t19 = g(3) * t43 + t39 * t32;
t40 = sin(qJ(3));
t44 = cos(qJ(3));
t49 = t39 * g(3) - t32 * t43;
t66 = t40 * t19 + t49 * t44;
t9 = t19 * t44 - t40 * t49;
t3 = t9 * t41 + t66 * t45;
t24 = -g(3) * pkin(4) + t32 * pkin(6);
t25 = t32 * pkin(4) + g(3) * pkin(6);
t14 = -t24 * t43 + t39 * t25;
t13 = pkin(3) * g(3) + t14;
t15 = t24 * t39 + t25 * t43;
t7 = pkin(3) * t32 + t15;
t68 = -t13 * t40 + t7 * t44;
t20 = g(3) * t44 + t32 * t40;
t67 = t20 * t45;
t61 = g(3) * t40;
t38 = sin(qJ(5));
t4 = -t41 * t66 + t9 * t45;
t59 = t4 * t38;
t57 = t13 * t44 + t40 * t7;
t31 = -g(1) * t56 + g(2) * t46;
t55 = t38 * t31;
t53 = t40 * t41;
t42 = cos(qJ(5));
t52 = t42 * t31;
t50 = g(3) * t45 + t32 * t41;
t48 = g(3) * pkin(2);
t35 = t44 * pkin(3) + pkin(2);
t33 = -pkin(4) * t39 + pkin(6) * t43;
t30 = pkin(4) * t43 + pkin(6) * t39 + pkin(3);
t29 = t40 * t45 + t44 * t41;
t28 = -t40 * t39 + t43 * t44;
t27 = t39 * t44 + t43 * t40;
t26 = -t44 * t45 + t53;
t18 = t32 * t44 - t61;
t16 = -t41 * t27 + t28 * t45;
t11 = t18 * t41 + t67;
t10 = t18 * t45 - t20 * t41;
t1 = t42 * t4;
t2 = [0, 0, 0, 0, 0, 0, -t31, t32, 0, 0, 0, 0, 0, 0, 0, 0, -t31 * t45, t31 * t41, t32, -pkin(1) * t31, 0, 0, 0, 0, 0, 0, t31 * t26, t31 * t29, t32, -t31 * (t45 * pkin(2) + pkin(1)), 0, 0, 0, 0, 0, 0, t31 * (t43 * t26 + t39 * t29), t31 * (-t26 * t39 + t29 * t43), t32, -t31 * (-pkin(3) * t53 + t35 * t45 + pkin(1)), 0, 0, 0, 0, 0, 0, -t16 * t52 + t38 * t32, t16 * t55 + t42 * t32, -t31 * (t27 * t45 + t41 * t28), -t31 * ((t30 * t44 + t33 * t40 + pkin(2)) * t45 + pkin(1) + (-t40 * t30 + t33 * t44) * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -g(3) * t41 + t32 * t45, 0, 0, 0, 0, 0, 0, 0, 0, t11, t10, 0, pkin(2) * t50, 0, 0, 0, 0, 0, 0, t4, -t3, 0, (pkin(2) * t51 + t35 * t62) * t41 + t45 * t48 + ((t44 * t51 - t61) * t41 + t67) * pkin(3), 0, 0, 0, 0, 0, 0, t1, -t59, t3, (pkin(2) * t32 + t68) * t41 + t45 * (t48 + t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t10, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, t11 * pkin(3), 0, 0, 0, 0, 0, 0, t1, -t59, t3, t41 * t68 + t57 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t59, t3, (t14 * t44 + t15 * t40) * t45 + (-t14 * t40 + t15 * t44) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t38 - t52, -t3 * t42 + t55, 0, 0;];
taug_reg = t2;
