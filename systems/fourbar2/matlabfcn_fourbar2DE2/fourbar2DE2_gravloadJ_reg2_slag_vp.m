% Calculate inertial parameters regressor of gravitation load for
% fourbar2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% 
% Output:
% taug_reg [1x(1*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:27
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbar2DE2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2DE2_gravloadJ_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2DE2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2DE2_gravloadJ_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t4 = cos(qJ(1));
t3 = sin(qJ(1));
t2 = g(1) * t4 + g(2) * t3;
t1 = g(1) * t3 - g(2) * t4;
t5 = [0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pkin(2) * t1, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t5;
