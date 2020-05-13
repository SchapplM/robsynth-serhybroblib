% Calculate potential energy for
% palh1m1DE2
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
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1DE2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m1DE2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_energypot_floatb_twist_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE2_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1DE2_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-14 20:02:04
% EndTime: 2020-04-14 20:02:11
% DurationCPUTime: 4.28s
% Computational Cost: add. (70670->173), mult. (105954->191), div. (5052->9), fcn. (67021->48), ass. (0->121)
t168 = -pkin(12) * m(6) + mrSges(5,2) - mrSges(6,3);
t81 = sin(qJ(4));
t87 = cos(qJ(4));
t167 = -pkin(10) * m(6) - mrSges(6,1) * t87 + mrSges(6,2) * t81 - mrSges(5,1);
t166 = -m(3) - m(9);
t165 = -m(4) - m(8);
t164 = -m(5) - m(6);
t163 = -m(2) - m(7) - m(10);
t161 = -m(11) - m(1) + t163;
t160 = mrSges(6,1) * t81 + mrSges(6,2) * t87 - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t109 = pkin(1) ^ 2;
t83 = sin(qJ(2));
t85 = sin(pkin(19));
t89 = cos(qJ(2));
t91 = cos(pkin(19));
t57 = t83 * t91 - t85 * t89;
t143 = pkin(7) * t57;
t133 = -0.2e1 * pkin(1) * t143 + t109;
t126 = pkin(7) ^ 2 + t133;
t47 = 0.1e1 / t126;
t136 = 0.1e1 / pkin(3) * t47;
t88 = cos(qJ(3));
t154 = -t88 / 0.2e1;
t149 = -pkin(8) + pkin(3);
t150 = -pkin(8) - pkin(3);
t42 = sqrt(-((pkin(7) - t149) * (pkin(7) + t149) + t133) * ((pkin(7) - t150) * (pkin(7) + t150) + t133));
t58 = t83 * t85 + t89 * t91;
t139 = t42 * t58;
t132 = -pkin(3) ^ 2 + pkin(8) ^ 2;
t45 = t126 - t132;
t50 = pkin(1) - t143;
t38 = -pkin(7) * t139 + t45 * t50;
t158 = -t38 / 0.2e1;
t39 = pkin(7) * t45 * t58 + t42 * t50;
t82 = sin(qJ(3));
t35 = (t154 * t39 + t158 * t82) * t136;
t157 = t39 / 0.2e1;
t36 = (t154 * t38 + t157 * t82) * t136;
t72 = pkin(23) + pkin(22);
t67 = sin(t72);
t68 = cos(t72);
t17 = t35 * t68 + t36 * t67;
t146 = pkin(5) * t17;
t134 = -0.2e1 * pkin(4) * t146 + pkin(5) ^ 2;
t124 = pkin(4) ^ 2 + t134;
t12 = 0.1e1 / t124;
t140 = t12 / pkin(11);
t155 = sin(pkin(21)) / 0.2e1;
t138 = pkin(9) ^ 2 - pkin(11) ^ 2;
t11 = t124 - t138;
t14 = -pkin(4) * t17 + pkin(5);
t18 = -t35 * t67 + t36 * t68;
t147 = -pkin(9) + pkin(11);
t148 = -pkin(9) - pkin(11);
t9 = sqrt(-((pkin(4) - t147) * (pkin(4) + t147) + t134) * ((pkin(4) - t148) * (pkin(4) + t148) + t134));
t142 = t18 * t9;
t7 = -pkin(4) * t142 + t11 * t14;
t74 = qJ(2) + qJ(3);
t79 = cos(pkin(21));
t8 = pkin(4) * t11 * t18 + t14 * t9;
t3 = atan2((t7 * t155 + t8 * t79 / 0.2e1) * t140, (-t7 * t79 / 0.2e1 + t8 * t155) * t140) + t74;
t1 = sin(t3);
t156 = sin(pkin(23)) / 0.2e1;
t78 = cos(pkin(23));
t22 = qJ(2) + atan2((t156 * t38 + t157 * t78) * t136, (t156 * t39 + t158 * t78) * t136);
t19 = pkin(22) - t22;
t2 = cos(t3);
t20 = sin(t22);
t21 = cos(t22);
t137 = 0.1e1 / pkin(8) * t47;
t153 = cos(pkin(18)) / 0.2e1;
t44 = t126 + t132;
t51 = pkin(1) * t57 - pkin(7);
t37 = -pkin(1) * t139 - t44 * t51;
t40 = pkin(1) * t44 * t58 - t42 * t51;
t86 = sin(pkin(18));
t26 = atan2((t40 * t153 + t37 * t86 / 0.2e1) * t137, (t37 * t153 - t86 * t40 / 0.2e1) * t137);
t23 = sin(t26);
t24 = cos(t26);
t107 = pkin(2) ^ 2;
t141 = 0.1e1 / pkin(2) / 0.2e1;
t122 = 0.1e1 / pkin(13) * t141;
t77 = sin(pkin(20));
t80 = cos(pkin(20));
t145 = pkin(6) * (-t77 * t88 - t80 * t82);
t130 = pkin(1) * t145;
t102 = pkin(6) ^ 2;
t131 = t102 + t109;
t49 = -0.2e1 * t130;
t125 = t49 + t131;
t123 = 0.1e1 / t125 * t141;
t144 = pkin(6) * (t77 * t82 - t80 * t88);
t135 = t102 + t49;
t151 = -pkin(2) + pkin(13);
t152 = -pkin(2) - pkin(13);
t41 = sqrt(-((pkin(1) - t151) * (pkin(1) + t151) + t135) * ((pkin(1) - t152) * (pkin(1) + t152) + t135));
t93 = pkin(13) ^ 2;
t43 = t107 - t93 + t125;
t48 = -pkin(1) + t145;
t33 = qJ(2) + atan2((t144 * t43 - t41 * t48) * t123, (-t144 * t41 - t43 * t48) * t123);
t29 = atan2(t41 * t122, (t107 + t93 + 0.2e1 * t130 - t131) * t122) + t33;
t27 = sin(t29);
t28 = cos(t29);
t31 = sin(t33);
t32 = cos(t33);
t10 = t124 + t138;
t127 = t12 / pkin(9) / 0.2e1;
t13 = -pkin(4) + t146;
t6 = -atan2((pkin(5) * t10 * t18 - t13 * t9) * t127, (-pkin(5) * t142 - t10 * t13) * t127) + t19;
t4 = sin(t6);
t5 = cos(t6);
t66 = -pkin(1) * t83 + pkin(16);
t69 = sin(t74);
t70 = cos(t74);
t159 = -mrSges(2,1) - m(11) * (pkin(4) * sin(t19) + t66) + t4 * mrSges(11,1) - t5 * mrSges(11,2) - m(10) * (-pkin(2) * t31 + pkin(16)) - t27 * mrSges(10,1) - t28 * mrSges(10,2) + m(7) * pkin(15) - t24 * mrSges(7,1) + t23 * mrSges(7,2) + mrSges(9,1) * t31 + mrSges(9,2) * t32 + mrSges(8,1) * t20 + mrSges(8,2) * t21 - mrSges(4,1) * t70 + mrSges(4,2) * t69 + mrSges(3,1) * t83 + mrSges(3,2) * t89 + t167 * t2 + t168 * t1;
t73 = pkin(14) + r_base(3);
t62 = t89 * pkin(1) + t73;
t90 = cos(qJ(1));
t84 = sin(qJ(1));
t59 = pkin(5) * t70 + t66;
t15 = (-mrSges(4,1) * t69 - mrSges(4,2) * t70 + m(7) * pkin(17) + mrSges(9,2) * t31 - t27 * mrSges(10,2) + t28 * mrSges(10,1) - t23 * mrSges(7,1) - t24 * mrSges(7,2) - mrSges(3,1) * t89 + mrSges(3,2) * t83 - mrSges(8,1) * t21 + mrSges(8,2) * t20 - m(11) * pkin(4) * cos(t19) + t5 * mrSges(11,1) + t4 * mrSges(11,2) - m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(11) + t165) * t62 + (-pkin(2) * m(10) - mrSges(9,1)) * t32 + t164 * (pkin(5) * t69 + t62) - t168 * t2 + t167 * t1 + (t163 + t166) * t73) * g(3) + (-mrSges(1,2) + t166 * (pkin(16) * t84 + r_base(2)) + t165 * (t66 * t84 + r_base(2)) + t164 * (t59 * t84 + r_base(2)) + t161 * r_base(2) + t160 * t90 + t159 * t84) * g(2) + (-mrSges(1,1) + t166 * (pkin(16) * t90 + r_base(1)) + t165 * (t66 * t90 + r_base(1)) + t164 * (t59 * t90 + r_base(1)) + t161 * r_base(1) + t159 * t90 - t160 * t84) * g(1);
U = t15;
