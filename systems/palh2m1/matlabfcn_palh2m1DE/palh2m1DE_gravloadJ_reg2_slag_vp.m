% Calculate inertial parameters regressor of gravitation load for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh2m1DE_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t21 = cos(qJ(2));
t18 = sin(qJ(1));
t28 = g(2) * t18;
t22 = cos(qJ(1));
t29 = g(1) * t22;
t11 = t28 + t29;
t16 = sin(qJ(3));
t20 = cos(qJ(3));
t8 = t20 * g(3) + t16 * t11;
t31 = t8 * t21;
t27 = t16 * g(3);
t25 = t21 * pkin(2);
t17 = sin(qJ(2));
t24 = t17 * t16;
t13 = t20 * pkin(3) + pkin(2);
t23 = -pkin(3) * t24 + t13 * t21 + pkin(1);
t19 = cos(qJ(4));
t15 = sin(qJ(4));
t10 = g(1) * t18 - g(2) * t22;
t9 = t11 * t20 - t27;
t7 = g(3) * t21 + t11 * t17;
t6 = t10 * t19 + t15 * t11;
t5 = -t15 * t10 + t19 * t11;
t4 = t9 * t17 + t31;
t3 = -t8 * t17 + t9 * t21;
t2 = t4 * pkin(3);
t1 = (pkin(2) * t28 + t13 * t29) * t17 + g(3) * t25 + ((t20 * t28 - t27) * t17 + t31) * pkin(3);
t12 = [0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t21, -t17 * t10, t11, pkin(1) * t10, 0, 0, 0, 0, 0, 0, t10 * (t21 * t20 - t24), -t10 * (t21 * t16 + t20 * t17), t11, t10 * (pkin(1) + t25), 0, 0, 0, 0, 0, 0, t10, 0, t11, t23 * t10, 0, 0, 0, 0, 0, 0, t6, t5, 0, t10 * (pkin(4) + t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -g(3) * t17 + t21 * t11, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, pkin(2) * t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0;];
taug_reg = t12;
