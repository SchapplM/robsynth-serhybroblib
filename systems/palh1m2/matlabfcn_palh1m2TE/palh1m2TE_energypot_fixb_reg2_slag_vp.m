% Calculate inertial parameters regressor of potential energy for
% palh1m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = palh1m2TE_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2TE_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_energypot_fixb_reg2_slag_vp: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:47:05
% EndTime: 2020-05-01 20:47:05
% DurationCPUTime: 0.51s
% Computational Cost: add. (331->89), mult. (754->148), div. (0->0), fcn. (803->20), ass. (0->74)
t153 = sin(pkin(21));
t156 = cos(pkin(21));
t157 = cos(pkin(20));
t196 = sin(pkin(20));
t134 = -t157 * t153 - t196 * t156;
t135 = -t196 * t153 + t156 * t157;
t152 = sin(pkin(22));
t155 = cos(pkin(22));
t117 = t134 * t155 - t135 * t152;
t163 = sin(pkin(18));
t169 = cos(pkin(18));
t188 = t134 * t152 + t135 * t155;
t184 = -t117 * t163 + t169 * t188;
t162 = sin(qJ(1));
t168 = cos(qJ(1));
t186 = g(1) * t168 + g(2) * t162;
t206 = t117 * t169 + t163 * t188;
t207 = t184 * g(3) - t186 * t206;
t201 = g(3) * pkin(13);
t161 = sin(qJ(2));
t200 = pkin(1) * t161;
t197 = g(3) * t206;
t160 = sin(qJ(3));
t149 = pkin(5) * t160 + pkin(1);
t193 = t149 * t161;
t166 = cos(qJ(3));
t192 = t161 * t166;
t167 = cos(qJ(2));
t191 = t166 * t167;
t190 = -pkin(1) * t167 - pkin(13);
t187 = -pkin(5) * t191 - pkin(15);
t154 = sin(pkin(19));
t158 = cos(pkin(19));
t137 = t154 * t166 + t158 * t160;
t138 = -t154 * t160 + t158 * t166;
t182 = t137 * t167 + t138 * t161;
t181 = t137 * t161 - t138 * t167;
t164 = sin(pkin(17));
t170 = cos(pkin(17));
t139 = t163 * t170 - t169 * t164;
t140 = t164 * t163 + t169 * t170;
t180 = -t139 * t167 + t140 * t161;
t179 = t139 * t161 + t140 * t167;
t178 = t152 * t169 - t155 * t163;
t177 = t152 * t163 + t155 * t169;
t176 = t160 * t167 + t192;
t175 = t160 * t161 - t191;
t172 = -pkin(5) * t192 - t149 * t167 - pkin(13);
t165 = cos(qJ(4));
t159 = sin(qJ(4));
t151 = g(1) * t200;
t150 = g(2) * t200;
t145 = -g(1) * t162 + g(2) * t168;
t144 = pkin(9) * t157 - pkin(11) * t196;
t143 = t196 * pkin(9) + pkin(11) * t157;
t142 = g(2) * t193;
t141 = g(1) * t193;
t136 = -t186 * pkin(15) - t201;
t133 = t138 * pkin(2);
t132 = t137 * pkin(2);
t131 = (-t152 * t153 + t155 * t156) * pkin(4);
t130 = (t152 * t156 + t153 * t155) * pkin(4);
t127 = g(3) * t161 + t186 * t167;
t126 = -g(3) * t167 + t186 * t161;
t125 = -t143 * t153 + t144 * t156;
t124 = t143 * t156 + t144 * t153;
t123 = (-pkin(15) * g(1) + t151) * t168 + (-g(2) * pkin(15) + t150) * t162 + t190 * g(3);
t116 = t175 * g(3) + t186 * t176;
t115 = -t176 * g(3) + t186 * t175;
t114 = -t124 * t152 + t125 * t155;
t113 = t124 * t155 + t125 * t152;
t111 = t184 * g(2);
t110 = t184 * g(1);
t1 = [0, 0, 0, 0, 0, 0, -t186, -t145, -g(3), -t201, 0, 0, 0, 0, 0, 0, t126, t127, t145, t136, 0, 0, 0, 0, 0, 0, t115, t116, t145, t123, 0, 0, 0, 0, 0, 0, t110 * t168 + t111 * t162 + t197, t207, t145, (t187 * g(1) + t141) * t168 + (t187 * g(2) + t142) * t162 + t172 * g(3), 0, 0, 0, 0, 0, 0, (g(2) * t159 + t110 * t165) * t168 + (-g(1) * t159 + t111 * t165) * t162 + t165 * t197, (g(2) * t165 - t110 * t159) * t168 + (-g(1) * t165 - t111 * t159) * t162 - t159 * t197, -t207, t141 * t168 + t142 * t162 + (-t113 * t169 + t114 * t163 + t172) * g(3) + t186 * (t113 * t163 + t114 * t169 + t187), 0, 0, 0, 0, 0, 0, -t179 * g(3) + t186 * t180, t180 * g(3) + t186 * t179, t145, -g(3) * (pkin(13) - pkin(16)) + t186 * pkin(14), 0, 0, 0, 0, 0, 0, -t178 * g(3) + t186 * t177, t177 * g(3) + t186 * t178, t145, t123, 0, 0, 0, 0, 0, 0, -t182 * g(3) + t186 * t181, t181 * g(3) + t186 * t182, t145, t136, 0, 0, 0, 0, 0, 0, t126, t127, t145, (-t132 * t167 - t133 * t161 - pkin(13)) * g(3) + t186 * (t132 * t161 - t133 * t167 - pkin(15)), 0, 0, 0, 0, 0, 0, t115, t116, t145, t150 * t162 + t151 * t168 + (-t130 * t169 + t131 * t163 + t190) * g(3) + t186 * (t130 * t163 + t131 * t169 - pkin(15));];
U_reg = t1;
