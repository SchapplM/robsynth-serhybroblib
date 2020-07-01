% Calculate minimal parameter regressor of potential energy for
% fourbar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% U_reg [1x9]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:43
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbar1OL_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1OL_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1OL_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1OL_energypot_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:43:18
% EndTime: 2020-06-26 17:43:18
% DurationCPUTime: 0.02s
% Computational Cost: add. (10->7), mult. (12->12), div. (0->0), fcn. (12->6), ass. (0->8)
t9 = cos(qJ(1));
t8 = cos(qJ(3));
t7 = sin(qJ(1));
t6 = sin(qJ(3));
t5 = qJ(1) + qJ(2);
t4 = cos(t5);
t3 = sin(t5);
t1 = [0, -g(1) * t9 - g(2) * t7, g(1) * t7 - g(2) * t9, 0, g(1) * t4 + g(2) * t3, -g(1) * t3 + g(2) * t4, 0, -g(1) * t8 - g(2) * t6, g(1) * t6 - g(2) * t8;];
U_reg = t1;
