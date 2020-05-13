% Calculate inertial parameters regressor of gravitation load for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh2m2DE_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t22 = sin(qJ(1));
t17 = cos(qJ(2));
t21 = t17 * pkin(4) + pkin(1);
t20 = pkin(2) + t21;
t16 = cos(qJ(3));
t19 = t16 * pkin(5) + t20;
t18 = cos(qJ(1));
t8 = g(1) * t18 + g(2) * t22;
t15 = cos(qJ(4));
t14 = sin(qJ(2));
t13 = sin(qJ(3));
t12 = sin(qJ(4));
t7 = -g(1) * t22 + g(2) * t18;
t6 = -t16 * g(3) + t13 * t8;
t5 = -t17 * g(3) + t14 * t8;
t4 = pkin(4) * t5;
t3 = pkin(5) * t6;
t2 = t12 * t8 - t15 * t7;
t1 = t12 * t7 + t15 * t8;
t9 = [0, 0, 0, 0, 0, 0, -t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t17, t14 * t7, -t8, -pkin(1) * t7, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, -t7 * t21, 0, 0, 0, 0, 0, 0, -t7 * t16, t13 * t7, -t8, -t20 * t7, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, -t19 * t7, 0, 0, 0, 0, 0, 0, t2, t1, 0, -(pkin(3) + t19) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, g(3) * t14 + t17 * t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, g(3) * t13 + t16 * t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t9;
