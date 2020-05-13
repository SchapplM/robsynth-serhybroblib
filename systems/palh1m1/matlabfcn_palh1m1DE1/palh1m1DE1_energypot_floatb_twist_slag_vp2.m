% Calculate potential energy for
% palh1m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1DE1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m1DE1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_energypot_floatb_twist_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE1_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1DE1_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 14:47:33
% EndTime: 2020-04-13 14:47:41
% DurationCPUTime: 7.69s
% Computational Cost: add. (159842->192), mult. (240243->252), div. (11448->9), fcn. (152251->46), ass. (0->140)
t105 = sin(qJ(2));
t111 = cos(qJ(2));
t124 = pkin(6) ^ 2;
t131 = pkin(1) ^ 2;
t149 = t124 + t131;
t160 = 0.1e1 / pkin(2) / 0.2e1;
t102 = cos(pkin(20));
t104 = sin(qJ(3));
t110 = cos(qJ(3));
t98 = sin(pkin(20));
t165 = pkin(6) * (-t102 * t104 - t110 * t98);
t147 = pkin(1) * t165;
t72 = -0.2e1 * t147;
t143 = 0.1e1 / (t72 + t149) * t160;
t164 = pkin(6) * (-t102 * t110 + t104 * t98);
t153 = t124 + t72;
t167 = pkin(13) - pkin(2);
t172 = -pkin(2) - pkin(13);
t58 = sqrt(-((pkin(1) - t167) * (pkin(1) + t167) + t153) * ((pkin(1) - t172) * (pkin(1) + t172) + t153));
t129 = pkin(2) ^ 2;
t142 = -pkin(13) ^ 2 + t149;
t60 = t129 + t72 + t142;
t71 = -pkin(1) + t165;
t48 = atan2((t164 * t60 - t58 * t71) * t143, (-t164 * t58 - t60 * t71) * t143);
t46 = sin(t48);
t47 = cos(t48);
t137 = t105 * t46 - t111 * t47;
t141 = 0.1e1 / pkin(13) * t160;
t57 = atan2(t58 * t141, (t129 - t142 + 0.2e1 * t147) * t141);
t55 = sin(t57);
t56 = cos(t57);
t180 = -m(10) * pkin(2) + t56 * mrSges(10,1) - t55 * mrSges(10,2) - mrSges(9,1);
t182 = t55 * mrSges(10,1) + t56 * mrSges(10,2) - mrSges(9,2);
t32 = -t105 * t47 - t111 * t46;
t108 = sin(pkin(18));
t107 = sin(pkin(19));
t113 = cos(pkin(19));
t81 = t105 * t113 - t107 * t111;
t163 = pkin(7) * t81;
t151 = -0.2e1 * pkin(1) * t163 + t131;
t146 = pkin(7) ^ 2 + t151;
t64 = 0.1e1 / t146;
t155 = 0.1e1 / pkin(8) * t64;
t161 = cos(pkin(18)) / 0.2e1;
t170 = -pkin(8) + pkin(3);
t171 = -pkin(3) - pkin(8);
t59 = sqrt(-((pkin(7) - t170) * (pkin(7) + t170) + t151) * ((pkin(7) - t171) * (pkin(7) + t171) + t151));
t84 = t105 * t107 + t111 * t113;
t157 = t59 * t84;
t148 = pkin(3) ^ 2 - pkin(8) ^ 2;
t61 = t146 - t148;
t74 = pkin(1) * t81 - pkin(7);
t51 = -pkin(1) * t157 - t61 * t74;
t54 = pkin(1) * t61 * t84 - t59 * t74;
t45 = atan2((t54 * t161 + t51 * t108 / 0.2e1) * t155, (t51 * t161 - t108 * t54 / 0.2e1) * t155);
t42 = sin(t45);
t43 = cos(t45);
t191 = m(7) * pkin(15) + mrSges(3,1) * t105 - t43 * mrSges(7,1) + mrSges(3,2) * t111 + t42 * mrSges(7,2) + t137 * t182 + t32 * t180 - mrSges(2,1);
t190 = -m(2) - m(7);
t189 = -m(6) - m(5);
t188 = pkin(4) * m(11);
t187 = -m(3) - m(10) - m(9);
t186 = -m(1) + t190;
t185 = -m(11) - m(4) - m(8);
t183 = -m(6) * pkin(12) + mrSges(5,2) - mrSges(6,3);
t103 = sin(qJ(4));
t109 = cos(qJ(4));
t179 = -m(6) * pkin(10) - t109 * mrSges(6,1) + t103 * mrSges(6,2) - mrSges(5,1);
t177 = t103 * mrSges(6,1) + t109 * mrSges(6,2) - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t62 = t146 + t148;
t73 = pkin(1) - t163;
t52 = -pkin(7) * t157 + t62 * t73;
t176 = -t52 / 0.2e1;
t53 = pkin(7) * t62 * t84 + t59 * t73;
t175 = t53 / 0.2e1;
t174 = sin(pkin(23)) / 0.2e1;
t173 = sin(pkin(21)) / 0.2e1;
t169 = -pkin(9) - pkin(11);
t168 = -pkin(9) + pkin(11);
t154 = 0.1e1 / pkin(3) * t64;
t162 = -t110 / 0.2e1;
t49 = (t104 * t176 + t53 * t162) * t154;
t50 = (t104 * t175 + t162 * t52) * t154;
t93 = pkin(23) + pkin(22);
t88 = sin(t93);
t89 = cos(t93);
t38 = t49 * t89 + t50 * t88;
t166 = pkin(5) * t38;
t159 = pkin(1) * t105;
t152 = -0.2e1 * pkin(4) * t166 + pkin(5) ^ 2;
t19 = sqrt(-((pkin(4) - t168) * (pkin(4) + t168) + t152) * ((pkin(4) - t169) * (pkin(4) + t169) + t152));
t39 = -t49 * t88 + t50 * t89;
t158 = t19 * t39;
t145 = pkin(4) ^ 2 + t152;
t34 = 0.1e1 / t145;
t156 = 0.1e1 / pkin(11) * t34;
t150 = pkin(9) ^ 2 - pkin(11) ^ 2;
t94 = pkin(14) + r_base(3);
t106 = sin(qJ(1));
t86 = t106 * pkin(16) + r_base(2);
t112 = cos(qJ(1));
t87 = t112 * pkin(16) + r_base(1);
t144 = 0.1e1 / pkin(9) * t34 / 0.2e1;
t85 = t111 * pkin(1) + t94;
t99 = cos(pkin(23));
t44 = atan2((t174 * t52 + t175 * t99) * t154, (t53 * t174 + t176 * t99) * t154);
t40 = sin(t44);
t41 = cos(t44);
t24 = -t105 * t41 - t111 * t40;
t138 = t105 * t40 - t111 * t41;
t82 = t104 * t111 + t105 * t110;
t83 = -t104 * t105 + t110 * t111;
t77 = -t106 * t159 + t86;
t78 = -t112 * t159 + t87;
t30 = t145 + t150;
t35 = -pkin(4) + t166;
t133 = atan2((pkin(5) * t30 * t39 - t19 * t35) * t144, (-pkin(5) * t158 - t30 * t35) * t144);
t132 = sin(t133);
t101 = cos(pkin(21));
t100 = cos(pkin(22));
t96 = sin(pkin(22));
t70 = t82 * t112;
t69 = t83 * t112;
t68 = t82 * t106;
t67 = t83 * t106;
t36 = -pkin(4) * t38 + pkin(5);
t31 = t145 - t150;
t23 = t24 * t112;
t22 = t138 * t112;
t21 = t24 * t106;
t20 = t138 * t106;
t18 = pkin(4) * t31 * t39 + t19 * t36;
t17 = -pkin(4) * t158 + t31 * t36;
t16 = cos(t133);
t14 = atan2((t17 * t173 + t18 * t101 / 0.2e1) * t156, (-t17 * t101 / 0.2e1 + t18 * t173) * t156);
t13 = cos(t14);
t12 = sin(t14);
t11 = -t100 * t16 - t132 * t96;
t10 = -t100 * t132 + t16 * t96;
t1 = (m(7) * pkin(17) - (-t100 * t138 - t24 * t96) * t188 + t105 * mrSges(3,2) - t43 * mrSges(7,2) - mrSges(3,1) * t111 - t42 * mrSges(7,1) - (t10 * t24 - t11 * t138) * mrSges(11,1) - (t10 * t138 + t11 * t24) * mrSges(11,2) - m(1) * r_base(3) - mrSges(4,1) * t82 - mrSges(4,2) * t83 + mrSges(8,1) * t138 - mrSges(8,2) * t24 - mrSges(1,3) - mrSges(2,3) + t185 * t85 + t179 * (t12 * t83 + t13 * t82) + t189 * (t82 * pkin(5) + t85) + t183 * (t12 * t82 - t83 * t13) + t182 * t32 - t180 * t137 + (t187 + t190) * t94) * g(3) + (-(t10 * t20 + t11 * t21) * mrSges(11,1) - (-t10 * t21 + t11 * t20) * mrSges(11,2) - (t100 * t21 - t20 * t96) * t188 - mrSges(8,1) * t21 - mrSges(8,2) * t20 - mrSges(4,1) * t67 + mrSges(4,2) * t68 - mrSges(1,2) + t187 * t86 + t185 * t77 + t189 * (t67 * pkin(5) + t77) + t183 * (t12 * t67 + t68 * t13) + t186 * r_base(2) + t179 * (-t12 * t68 + t13 * t67)) * g(2) + (-(t100 * t23 - t22 * t96) * t188 - (t10 * t22 + t11 * t23) * mrSges(11,1) - (-t10 * t23 + t11 * t22) * mrSges(11,2) - mrSges(4,1) * t69 + mrSges(4,2) * t70 - mrSges(8,1) * t23 - mrSges(8,2) * t22 - mrSges(1,1) + t187 * t87 + t185 * t78 + t189 * (t69 * pkin(5) + t78) + t183 * (t12 * t69 + t70 * t13) + t186 * r_base(1) + t179 * (-t12 * t70 + t13 * t69)) * g(1) + (g(1) * t191 + t177 * g(2)) * t112 + (-t177 * g(1) + t191 * g(2)) * t106;
U = t1;
