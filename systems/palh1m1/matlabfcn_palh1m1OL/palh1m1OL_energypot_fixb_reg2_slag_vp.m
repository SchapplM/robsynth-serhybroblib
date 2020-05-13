% Calculate inertial parameters regressor of potential energy for
% palh1m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% 
% Output:
% U_reg [1x(13*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh1m1OL_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_energypot_fixb_reg2_slag_vp: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_energypot_fixb_reg2_slag_vp: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:36:40
% EndTime: 2020-04-15 19:36:40
% DurationCPUTime: 0.22s
% Computational Cost: add. (211->64), mult. (208->76), div. (0->0), fcn. (184->22), ass. (0->43)
t125 = sin(qJ(1));
t129 = cos(qJ(1));
t143 = -g(1) * t129 - g(2) * t125;
t142 = g(3) * pkin(13);
t121 = qJ(2) + qJ(3);
t116 = qJ(4) + t121;
t105 = sin(t116);
t139 = g(3) * t105;
t128 = cos(qJ(2));
t138 = t128 * pkin(1) + pkin(13);
t123 = sin(qJ(5));
t137 = t125 * t123;
t127 = cos(qJ(5));
t136 = t125 * t127;
t135 = t129 * t123;
t134 = t129 * t127;
t120 = qJ(2) + qJ(7);
t119 = qJ(2) + qJ(8);
t111 = sin(t121);
t133 = pkin(5) * t111 + t138;
t124 = sin(qJ(2));
t132 = -t124 * pkin(1) + pkin(15);
t108 = pkin(19) - t120;
t126 = cos(qJ(6));
t122 = sin(qJ(6));
t115 = qJ(9) + t119;
t114 = cos(t121);
t113 = cos(t120);
t112 = cos(t119);
t110 = sin(t120);
t109 = sin(t119);
t107 = cos(t116);
t106 = cos(t115);
t104 = sin(t115);
t101 = -qJ(10) + t108;
t100 = cos(t101);
t99 = sin(t101);
t97 = g(1) * t125 - g(2) * t129;
t96 = pkin(5) * t114 + t132;
t94 = pkin(15) * t143 - t142;
t93 = -g(3) * t107 - t143 * t105;
t92 = -g(3) * t138 + t132 * t143;
t1 = [0, 0, 0, 0, 0, 0, t143, t97, -g(3), -t142, 0, 0, 0, 0, 0, 0, -g(3) * t128 - t124 * t143, g(3) * t124 - t128 * t143, -t97, t94, 0, 0, 0, 0, 0, 0, -g(3) * t111 + t114 * t143, -g(3) * t114 - t111 * t143, -t97, t92, 0, 0, 0, 0, 0, 0, t143 * t107 - t139, t93, -t97, -g(3) * t133 + t143 * t96, 0, 0, 0, 0, 0, 0, -g(1) * (t107 * t134 + t137) - g(2) * (t107 * t136 - t135) - t127 * t139, -g(1) * (-t107 * t135 + t136) - g(2) * (-t107 * t137 - t134) + t123 * t139, -t93, -g(3) * (t105 * pkin(9) - t107 * pkin(11) + t133) + t143 * (pkin(9) * t107 + pkin(11) * t105 + t96), 0, 0, 0, 0, 0, 0, -g(3) * t122 + t126 * t143, -g(3) * t126 - t122 * t143, -t97, -g(3) * (-pkin(16) + pkin(13)) - t143 * pkin(14), 0, 0, 0, 0, 0, 0, -g(3) * t113 - t110 * t143, g(3) * t110 - t113 * t143, -t97, t92, 0, 0, 0, 0, 0, 0, -g(3) * t112 - t109 * t143, g(3) * t109 - t112 * t143, -t97, t94, 0, 0, 0, 0, 0, 0, g(3) * t106 + t104 * t143, -g(3) * t104 + t106 * t143, -t97, -g(3) * (pkin(2) * t112 + pkin(13)) + t143 * (-pkin(2) * t109 + pkin(15)), 0, 0, 0, 0, 0, 0, g(3) * t100 - t143 * t99, g(3) * t99 + t100 * t143, -t97, -g(3) * (pkin(4) * cos(t108) + t138) + t143 * (pkin(4) * sin(t108) + t132);];
U_reg = t1;
