% Calculate inertial parameters regressor of potential energy for
% palh3m2OL
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
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh3m2OL_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_energypot_fixb_reg2_slag_vp: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2OL_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_energypot_fixb_reg2_slag_vp: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:40:29
% EndTime: 2020-05-07 04:40:29
% DurationCPUTime: 0.16s
% Computational Cost: add. (173->54), mult. (171->64), div. (0->0), fcn. (153->18), ass. (0->36)
t103 = cos(qJ(1));
t99 = sin(qJ(1));
t105 = g(1) * t103 + g(2) * t99;
t115 = g(3) * pkin(11);
t95 = qJ(2) + qJ(3);
t90 = qJ(4) + t95;
t82 = sin(t90);
t114 = g(3) * t82;
t97 = sin(qJ(5));
t112 = t99 * t97;
t98 = sin(qJ(2));
t111 = t98 * pkin(1) + pkin(11);
t102 = cos(qJ(2));
t110 = t102 * pkin(1) + pkin(12);
t109 = t103 * t97;
t101 = cos(qJ(5));
t108 = t99 * t101;
t107 = t103 * t101;
t94 = qJ(2) + qJ(7);
t89 = pkin(15) - t94;
t86 = sin(t95);
t106 = -pkin(4) * t86 + t111;
t100 = cos(qJ(6));
t96 = sin(qJ(6));
t88 = cos(t95);
t87 = cos(t94);
t85 = sin(t94);
t84 = -qJ(8) + t89;
t83 = cos(t90);
t80 = cos(t84);
t79 = sin(t84);
t78 = g(1) * t99 - g(2) * t103;
t77 = -pkin(4) * t88 + t110;
t75 = -g(3) * t83 + t105 * t82;
t74 = -g(3) * t111 - t105 * t110;
t1 = [0, 0, 0, 0, 0, 0, -t105, t78, -g(3), -t115, 0, 0, 0, 0, 0, 0, -g(3) * t98 - t105 * t102, -g(3) * t102 + t105 * t98, -t78, -t105 * pkin(12) - t115, 0, 0, 0, 0, 0, 0, g(3) * t86 + t105 * t88, g(3) * t88 - t105 * t86, -t78, t74, 0, 0, 0, 0, 0, 0, t105 * t83 + t114, -t75, -t78, -g(3) * t106 - t105 * t77, 0, 0, 0, 0, 0, 0, -g(1) * (-t83 * t107 + t112) - g(2) * (-t83 * t108 - t109) + t101 * t114, -g(1) * (t83 * t109 + t108) - g(2) * (t83 * t112 - t107) - t97 * t114, t75, -g(3) * (-t82 * pkin(8) + t83 * pkin(10) + t106) + t105 * (pkin(8) * t83 + pkin(10) * t82 - t77), 0, 0, 0, 0, 0, 0, -g(3) * t96 - t105 * t100, -g(3) * t100 + t105 * t96, -t78, -g(3) * (pkin(13) + pkin(11)) + t105 * pkin(6), 0, 0, 0, 0, 0, 0, -g(3) * t85 - t105 * t87, -g(3) * t87 + t105 * t85, -t78, t74, 0, 0, 0, 0, 0, 0, -g(3) * t79 + t105 * t80, g(3) * t80 + t105 * t79, -t78, -g(3) * (-pkin(3) * sin(t89) + t111) - t105 * (pkin(3) * cos(t89) + t110);];
U_reg = t1;
