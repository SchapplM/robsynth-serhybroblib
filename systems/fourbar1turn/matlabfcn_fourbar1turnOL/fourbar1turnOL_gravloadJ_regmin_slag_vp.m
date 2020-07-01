% Calculate minimal parameter regressor of gravitation load for
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
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:56
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbar1turnOL_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:56:29
% EndTime: 2020-06-27 16:56:29
% DurationCPUTime: 0.07s
% Computational Cost: add. (40->12), mult. (68->22), div. (0->0), fcn. (68->8), ass. (0->14)
t11 = cos(qJ(1));
t8 = sin(qJ(1));
t13 = g(1) * t8 - g(2) * t11;
t12 = g(1) * t11 + g(2) * t8;
t10 = cos(qJ(2));
t9 = cos(qJ(4));
t7 = sin(qJ(2));
t6 = sin(qJ(4));
t5 = qJ(2) + qJ(3);
t4 = cos(t5);
t3 = sin(t5);
t2 = -g(3) * t3 - t12 * t4;
t1 = g(3) * t4 - t12 * t3;
t14 = [0, t13, t12, 0, 0, 0, 0, 0, t13 * t10, -t13 * t7, 0, 0, 0, 0, 0, -t13 * t4, t13 * t3, 0, 0, 0, 0, 0, t13 * t9, -t13 * t6; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t10 + t12 * t7, g(3) * t7 + t12 * t10, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t9 + t12 * t6, g(3) * t6 + t12 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t14;
