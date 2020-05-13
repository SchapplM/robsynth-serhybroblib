% Calculate inertial parameters regressor of potential energy for
% fourbar2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% 
% Output:
% U_reg [1x(1*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:17
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbar2TE_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2TE_energypot_fixb_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2TE_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2TE_energypot_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:17:26
% EndTime: 2020-04-24 20:17:26
% DurationCPUTime: 0.02s
% Computational Cost: add. (12->11), mult. (12->6), div. (0->0), fcn. (10->2), ass. (0->5)
t6 = cos(qJ(1));
t5 = sin(qJ(1));
t4 = g(1) * t6 + g(2) * t5;
t3 = g(1) * t5 - g(2) * t6;
t1 = [0, 0, 0, 0, 0, 0, -t4, t3, -g(3), 0, 0, 0, 0, 0, 0, 0, -g(1), -g(2), -g(3), -pkin(2) * t4, 0, 0, 0, 0, 0, 0, -t4, t3, -g(3), -pkin(1) * g(1);];
U_reg = t1;
