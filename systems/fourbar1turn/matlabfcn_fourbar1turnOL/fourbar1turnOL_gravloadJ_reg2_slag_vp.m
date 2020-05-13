% Calculate inertial parameters regressor of gravitation load for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbar1turnOL_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t12 = cos(qJ(1));
t9 = sin(qJ(1));
t15 = g(1) * t9 - g(2) * t12;
t3 = g(1) * t12 + g(2) * t9;
t11 = cos(qJ(2));
t14 = t15 * t11;
t8 = sin(qJ(2));
t13 = -g(3) * t11 + t3 * t8;
t10 = cos(qJ(4));
t7 = sin(qJ(4));
t6 = qJ(2) + qJ(3);
t5 = cos(t6);
t4 = sin(t6);
t2 = -g(3) * t4 - t3 * t5;
t1 = g(3) * t5 - t3 * t4;
t16 = [0, 0, 0, 0, 0, 0, t15, t3, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15 * t8, -t3, 0, 0, 0, 0, 0, 0, 0, -t15 * t5, t15 * t4, -t3, pkin(2) * t14, 0, 0, 0, 0, 0, 0, t15 * t10, -t15 * t7, -t3, t15 * pkin(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, g(3) * t8 + t3 * t11, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t13 * pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t10 + t3 * t7, g(3) * t7 + t3 * t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t16;
