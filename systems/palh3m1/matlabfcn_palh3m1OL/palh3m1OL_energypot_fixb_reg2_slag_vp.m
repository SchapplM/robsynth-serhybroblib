% Calculate inertial parameters regressor of potential energy for
% palh3m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% U_reg [1x(10*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh3m1OL_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1OL_energypot_fixb_reg2_slag_vp: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1OL_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_energypot_fixb_reg2_slag_vp: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:10:55
% EndTime: 2020-04-20 17:10:55
% DurationCPUTime: 0.17s
% Computational Cost: add. (173->49), mult. (171->64), div. (0->0), fcn. (153->18), ass. (0->36)
t102 = cos(qJ(1));
t98 = sin(qJ(1));
t104 = g(1) * t102 + g(2) * t98;
t114 = g(3) * pkin(11);
t94 = qJ(2) + qJ(3);
t90 = qJ(4) + t94;
t82 = sin(t90);
t113 = g(3) * t82;
t96 = sin(qJ(5));
t111 = t98 * t96;
t97 = sin(qJ(2));
t110 = t97 * pkin(1) + pkin(11);
t101 = cos(qJ(2));
t109 = t101 * pkin(1) + pkin(12);
t108 = t102 * t96;
t100 = cos(qJ(5));
t107 = t98 * t100;
t106 = t102 * t100;
t93 = qJ(2) + qJ(7);
t89 = pkin(15) - t93;
t86 = sin(t94);
t105 = -pkin(4) * t86 + t110;
t99 = cos(qJ(6));
t95 = sin(qJ(6));
t88 = cos(t94);
t87 = cos(t93);
t85 = sin(t93);
t84 = -qJ(8) + t89;
t83 = cos(t90);
t80 = cos(t84);
t79 = sin(t84);
t78 = -g(1) * t98 + g(2) * t102;
t77 = -pkin(4) * t88 + t109;
t75 = g(3) * t83 - t104 * t82;
t74 = -g(3) * t110 - t104 * t109;
t1 = [0, 0, 0, 0, 0, 0, -t104, -t78, -g(3), -t114, 0, 0, 0, 0, 0, 0, -g(3) * t97 - t104 * t101, -g(3) * t101 + t104 * t97, t78, -t104 * pkin(12) - t114, 0, 0, 0, 0, 0, 0, g(3) * t86 + t104 * t88, g(3) * t88 - t104 * t86, t78, t74, 0, 0, 0, 0, 0, 0, t104 * t83 + t113, t75, t78, -g(3) * t105 - t104 * t77, 0, 0, 0, 0, 0, 0, t100 * t113 - g(2) * (-t83 * t107 - t108) - g(1) * (-t83 * t106 + t111), -t96 * t113 - g(2) * (t83 * t111 - t106) - g(1) * (t83 * t108 + t107), -t75, -g(3) * (-t82 * pkin(8) + t83 * pkin(10) + t105) + t104 * (pkin(8) * t83 + pkin(10) * t82 - t77), 0, 0, 0, 0, 0, 0, -g(3) * t95 - t104 * t99, -g(3) * t99 + t104 * t95, t78, -g(3) * (pkin(13) + pkin(11)) + t104 * pkin(6), 0, 0, 0, 0, 0, 0, -g(3) * t85 - t104 * t87, -g(3) * t87 + t104 * t85, t78, t74, 0, 0, 0, 0, 0, 0, -g(3) * t79 + t104 * t80, g(3) * t80 + t104 * t79, t78, -g(3) * (-pkin(3) * sin(t89) + t110) - t104 * (pkin(3) * cos(t89) + t109);];
U_reg = t1;
