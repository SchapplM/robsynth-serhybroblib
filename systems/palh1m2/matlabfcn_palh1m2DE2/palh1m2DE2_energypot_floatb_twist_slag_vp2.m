% Calculate potential energy for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
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
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m2DE2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m2DE2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_energypot_floatb_twist_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2DE2_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:57:49
% EndTime: 2020-05-02 20:57:51
% DurationCPUTime: 0.81s
% Computational Cost: add. (392->113), mult. (284->110), div. (0->0), fcn. (219->45), ass. (0->59)
t39 = cos(qJ(2));
t70 = g(3) * t39;
t65 = m(5) + m(6);
t19 = t65 * pkin(5) + mrSges(4,1);
t32 = sin(qJ(3));
t36 = sin(pkin(17));
t38 = cos(qJ(3));
t42 = cos(pkin(17));
t35 = sin(pkin(18));
t41 = cos(pkin(18));
t46 = mrSges(7,1) * t35 - mrSges(7,2) * t41;
t47 = mrSges(7,1) * t41 + mrSges(7,2) * t35;
t67 = m(8) + m(4) + t65;
t2 = t67 * pkin(1) + mrSges(4,2) * t38 + t19 * t32 + t46 * t36 + t47 * t42 + mrSges(3,1);
t3 = -t32 * mrSges(4,2) + t19 * t38 - t47 * t36 + t46 * t42 - mrSges(3,2);
t33 = sin(qJ(2));
t55 = -m(9) - m(10) - m(3) - t67;
t69 = pkin(14) * m(7) + t2 * t33 - t3 * t39 - mrSges(2,1) + (-m(11) + t55) * pkin(15);
t68 = g(1) * r_base(1) + r_base(2) * g(2);
t31 = sin(qJ(4));
t37 = cos(qJ(4));
t66 = mrSges(6,1) * t31 + mrSges(6,2) * t37 - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t63 = m(11) * g(1);
t62 = m(11) * g(2);
t60 = mrSges(11,1) * g(2);
t59 = mrSges(11,2) * g(1);
t58 = mrSges(11,2) * g(2);
t28 = sin(pkin(19));
t29 = cos(pkin(19));
t13 = t28 * t38 + t29 * t32;
t14 = t28 * t32 - t29 * t38;
t7 = qJ(2) + atan2(t14, t13);
t25 = pkin(18) - pkin(22);
t54 = -t63 / 0.2e1;
t53 = -t62 / 0.2e1;
t52 = t62 / 0.2e1;
t23 = -qJ(2) + t25;
t51 = pkin(21) - atan2(cos(t23), -sin(t23)) - qJ(2);
t49 = -pkin(21) + t25;
t11 = -qJ(1) + t51;
t10 = qJ(1) + t51;
t34 = sin(qJ(1));
t40 = cos(qJ(1));
t48 = g(1) * t40 + g(2) * t34;
t44 = m(7) + m(2) - t55;
t30 = mrSges(11,1) * g(1);
t27 = qJ(1) - qJ(2);
t26 = qJ(1) + qJ(2);
t24 = pkin(2) * m(10) + mrSges(9,1);
t21 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t20 = -pkin(20) + t49;
t18 = -qJ(2) - qJ(3) + t49;
t15 = pkin(9) * m(6) + mrSges(6,1) * t37 - mrSges(6,2) * t31 + mrSges(5,1);
t9 = atan2(-sin(t18), cos(t18));
t6 = -t9 + t51;
t5 = -t9 + t11;
t4 = -t9 + t10;
t1 = atan2(t14, -t13) + t7;
t8 = (t30 - t58) * sin(t4) / 0.2e1 + (t30 + t58) * sin(t5) / 0.2e1 + (-t59 - t60) * cos(t4) / 0.2e1 + (-t59 + t60) * cos(t5) / 0.2e1 + (t48 * mrSges(8,1) + g(3) * mrSges(8,2)) * cos(t25) + (-g(3) * t21 + t48 * t15) * cos(t20) + (g(3) * t15 + t48 * t21) * sin(t20) + (mrSges(9,2) * g(3) + t48 * t24) * sin(t7) + (t48 * mrSges(9,2) - g(3) * t24) * cos(t7) + (g(3) * mrSges(8,1) - t48 * mrSges(8,2)) * sin(t25) + (-t48 * mrSges(10,1) - g(3) * mrSges(10,2)) * sin(t1) + (g(3) * mrSges(10,1) - t48 * mrSges(10,2)) * cos(t1) - (-pkin(16) * m(7) + t44 * pkin(13) + mrSges(1,3) + mrSges(2,3)) * g(3) - g(3) * t3 * t33 + mrSges(11,1) * cos(t6) * g(3) + mrSges(11,2) * sin(t6) * g(3) + ((sin(t11) + sin(t10)) * t54 - m(11) * g(3) * cos(t51) + cos(t10) * t52 + cos(t11) * t53) * pkin(4) + (-m(11) * t70 + cos(t26) * t53 + sin(t27) * t54 + cos(t27) * t52 + sin(t26) * t63 / 0.2e1) * pkin(1) - t2 * t70 - (g(3) * (pkin(13) + r_base(3)) + t68) * m(11) + (-g(3) * r_base(3) - t68) * (m(1) + t44) + (t69 * g(1) + t66 * g(2)) * t40 + (-t66 * g(1) + t69 * g(2)) * t34 - mrSges(1,1) * g(1) - mrSges(1,2) * g(2);
U = t8;
