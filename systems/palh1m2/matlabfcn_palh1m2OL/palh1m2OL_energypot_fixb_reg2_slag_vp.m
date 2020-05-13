% Calculate inertial parameters regressor of potential energy for
% palh1m2OL
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
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh1m2OL_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_energypot_fixb_reg2_slag_vp: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2OL_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_energypot_fixb_reg2_slag_vp: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:46:52
% EndTime: 2020-05-02 21:46:52
% DurationCPUTime: 0.34s
% Computational Cost: add. (226->76), mult. (389->110), div. (0->0), fcn. (328->24), ass. (0->54)
t137 = sin(qJ(2));
t145 = cos(qJ(2));
t146 = cos(qJ(1));
t154 = t146 * g(1);
t138 = sin(qJ(1));
t156 = t138 * g(2);
t116 = t154 + t156;
t135 = sin(qJ(4));
t143 = cos(qJ(4));
t106 = t143 * g(3) - t135 * t116;
t157 = t135 * g(3);
t111 = t116 * t143 + t157;
t136 = sin(qJ(3));
t144 = cos(qJ(3));
t97 = t106 * t144 - t111 * t136;
t99 = t106 * t136 + t111 * t144;
t162 = t137 * t97 + t99 * t145;
t161 = -t99 * t137 + t97 * t145;
t159 = pkin(13) * g(3);
t158 = pkin(5) * t144;
t155 = t144 * g(3);
t126 = pkin(11) * t135 + pkin(5);
t130 = qJ(2) + qJ(8);
t129 = pkin(19) - qJ(7);
t150 = pkin(9) * t156;
t149 = t136 * g(3) + t116 * t144;
t108 = -g(3) * t145 + t116 * t137;
t113 = -t116 * pkin(15) - t159;
t147 = pkin(9) * t157 + (pkin(9) * t143 + t126) * t154 + (-pkin(11) * g(3) + t150) * t143 + t126 * t156;
t142 = cos(qJ(5));
t141 = cos(qJ(6));
t140 = cos(qJ(7));
t139 = cos(qJ(8));
t134 = sin(qJ(5));
t133 = sin(qJ(6));
t132 = sin(qJ(7));
t131 = sin(qJ(8));
t127 = qJ(9) + t130;
t125 = t136 * pkin(5) + pkin(1);
t124 = cos(t127);
t123 = sin(t127);
t121 = -qJ(2) - qJ(10) + t129;
t120 = cos(t121);
t119 = sin(t121);
t117 = -g(1) * t138 + g(2) * t146;
t112 = g(3) * t137 + t116 * t145;
t110 = t132 * g(3) + t116 * t140;
t109 = t131 * g(3) + t116 * t139;
t107 = -t136 * t116 + t155;
t105 = t140 * g(3) - t132 * t116;
t104 = t139 * g(3) - t131 * t116;
t101 = -g(3) * (t145 * pkin(1) + pkin(13)) + t116 * (t137 * pkin(1) - pkin(15));
t100 = -(-pkin(9) * t135 + pkin(11) * t143) * t154 + (-pkin(9) * g(3) - pkin(11) * t156) * t143 + t135 * t150 - g(3) * t126;
t1 = [0, 0, 0, 0, 0, 0, -t116, -t117, -g(3), -t159, 0, 0, 0, 0, 0, 0, t108, t112, t117, t113, 0, 0, 0, 0, 0, 0, -t107 * t137 - t149 * t145, -t107 * t145 + t137 * t149, t117, t101, 0, 0, 0, 0, 0, 0, -t162, -t161, t117, (-pkin(5) * t155 + t116 * t125) * t137 - (t145 * t158 + pkin(15)) * t154 + (-t125 * g(3) - t156 * t158) * t145 - pkin(15) * t156 - t159, 0, 0, 0, 0, 0, 0, t134 * t117 - t162 * t142, t142 * t117 + t162 * t134, t161, (pkin(1) * t116 + t100 * t144 + t147 * t136) * t137 + (-pkin(1) * g(3) + t100 * t136 - t147 * t144) * t145 + t113, 0, 0, 0, 0, 0, 0, -t133 * g(3) - t116 * t141, -t141 * g(3) + t116 * t133, t117, -(pkin(13) - pkin(16)) * g(3) + t116 * pkin(14), 0, 0, 0, 0, 0, 0, -t105 * t145 + t137 * t110, t105 * t137 + t110 * t145, t117, t101, 0, 0, 0, 0, 0, 0, -t104 * t145 + t137 * t109, t104 * t137 + t109 * t145, t117, t113, 0, 0, 0, 0, 0, 0, g(3) * t124 - t116 * t123, -g(3) * t123 - t116 * t124, t117, (sin(t130) * t116 - cos(t130) * g(3)) * pkin(2) + t113, 0, 0, 0, 0, 0, 0, g(3) * t120 + t116 * t119, g(3) * t119 - t116 * t120, t117, (cos(t129) * t108 - sin(t129) * t112) * pkin(4) + t108 * pkin(1) + t113;];
U_reg = t1;
