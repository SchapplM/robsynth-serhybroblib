% Calculate inertial parameters regressor of potential energy for
% fourbar2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbar2OL_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2OL_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_energypot_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:26
% EndTime: 2020-04-24 20:32:26
% DurationCPUTime: 0.02s
% Computational Cost: add. (16->12), mult. (24->14), div. (0->0), fcn. (22->6), ass. (0->9)
t12 = cos(qJ(1));
t11 = cos(qJ(2));
t10 = cos(qJ(3));
t9 = sin(qJ(1));
t8 = sin(qJ(2));
t7 = sin(qJ(3));
t6 = g(1) * t12 + g(2) * t9;
t5 = g(1) * t9 - g(2) * t12;
t1 = [0, 0, 0, 0, 0, 0, -t6, t5, -g(3), 0, 0, 0, 0, 0, 0, 0, t6 * t11 - t5 * t8, -t11 * t5 - t6 * t8, -g(3), -pkin(2) * t6, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t7, g(1) * t7 - g(2) * t10, -g(3), -pkin(1) * g(1);];
U_reg = t1;
