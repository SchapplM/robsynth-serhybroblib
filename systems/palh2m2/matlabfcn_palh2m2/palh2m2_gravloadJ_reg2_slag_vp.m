% Calculate inertial parameters regressor of gravitation load for
% palh2m2
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
% Datum: 2019-06-06 14:46
% Revision: 7254ec7b167830f9592b38d39d95d449e6fd98ef (2019-06-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh2m2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t16 = sin(qJ(1));
t20 = cos(qJ(1));
t8 = g(1) * t16 - g(2) * t20;
t19 = cos(qJ(2));
t24 = t19 * pkin(4) + pkin(1);
t23 = pkin(2) + t24;
t9 = g(1) * t20 + g(2) * t16;
t14 = sin(qJ(3));
t18 = cos(qJ(3));
t22 = -g(3) * t18 + t9 * t14;
t15 = sin(qJ(2));
t21 = -g(3) * t19 + t9 * t15;
t17 = cos(qJ(4));
t13 = sin(qJ(4));
t7 = t18 * pkin(5) + t23;
t6 = -t20 * t13 - t16 * t17;
t5 = t16 * t13 - t20 * t17;
t4 = t21 * pkin(4);
t3 = t22 * pkin(5);
t2 = -g(1) * t6 + g(2) * t5;
t1 = -g(1) * t5 - g(2) * t6;
t10 = [0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t19, -t8 * t15, -t9, t8 * pkin(1), 0, 0, 0, 0, 0, 0, t8, 0, -t9, t8 * t24, 0, 0, 0, 0, 0, 0, t8 * t18, -t8 * t14, -t9, t8 * t23, 0, 0, 0, 0, 0, 0, t8, 0, -t9, t8 * t7, 0, 0, 0, 0, 0, 0, t2, t1, 0, t8 * (pkin(3) + t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, g(3) * t15 + t9 * t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, g(3) * t14 + t9 * t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg  = t10;
