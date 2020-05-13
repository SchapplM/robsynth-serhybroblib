% Calculate potential energy for
% palh1m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1OL_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_energypot_floatb_twist_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m1OL_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_energypot_floatb_twist_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:17
% EndTime: 2020-04-15 19:28:18
% DurationCPUTime: 0.70s
% Computational Cost: add. (284->99), mult. (246->76), div. (0->0), fcn. (184->22), ass. (0->41)
t65 = -pkin(11) * m(6) + mrSges(5,2) - mrSges(6,3);
t35 = sin(qJ(5));
t39 = cos(qJ(5));
t64 = -pkin(9) * m(6) - mrSges(6,1) * t39 + mrSges(6,2) * t35 - mrSges(5,1);
t63 = -m(3) - m(9);
t62 = -m(4) - m(8);
t61 = -m(5) - m(6);
t60 = -m(2) - m(7) - m(10);
t58 = -m(1) - m(11) + t60;
t32 = qJ(2) + qJ(7);
t20 = pkin(19) - t32;
t13 = -qJ(10) + t20;
t10 = cos(t13);
t36 = sin(qJ(2));
t15 = -t36 * pkin(1) + pkin(15);
t31 = qJ(2) + qJ(8);
t27 = qJ(9) + t31;
t16 = sin(t27);
t33 = qJ(2) + qJ(3);
t28 = qJ(4) + t33;
t17 = sin(t28);
t18 = cos(t27);
t19 = cos(t28);
t21 = sin(t31);
t22 = sin(t32);
t23 = sin(t33);
t24 = cos(t31);
t25 = cos(t32);
t26 = cos(t33);
t34 = sin(qJ(6));
t38 = cos(qJ(6));
t40 = cos(qJ(2));
t9 = sin(t13);
t57 = -mrSges(2,1) - m(11) * (pkin(4) * sin(t20) + t15) + t9 * mrSges(11,1) - t10 * mrSges(11,2) + m(7) * pkin(14) - t38 * mrSges(7,1) + t34 * mrSges(7,2) - m(10) * (-pkin(2) * t21 + pkin(15)) - t16 * mrSges(10,1) - t18 * mrSges(10,2) + mrSges(9,1) * t21 + mrSges(9,2) * t24 + mrSges(8,1) * t22 + mrSges(8,2) * t25 - mrSges(4,1) * t26 + mrSges(4,2) * t23 + mrSges(3,1) * t36 + mrSges(3,2) * t40 + t64 * t19 + t65 * t17;
t56 = t35 * mrSges(6,1) + t39 * mrSges(6,2) - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t30 = pkin(13) + r_base(3);
t7 = t40 * pkin(1) + t30;
t41 = cos(qJ(1));
t37 = sin(qJ(1));
t4 = pkin(5) * t26 + t15;
t1 = (-m(11) * pkin(4) * cos(t20) + m(7) * pkin(16) - mrSges(1,3) - mrSges(2,3) - m(1) * r_base(3) + t9 * mrSges(11,2) + t10 * mrSges(11,1) - t16 * mrSges(10,2) + t18 * mrSges(10,1) + t21 * mrSges(9,2) + t22 * mrSges(8,2) - t23 * mrSges(4,1) - t25 * mrSges(8,1) - t26 * mrSges(4,2) - t34 * mrSges(7,1) + t36 * mrSges(3,2) - t38 * mrSges(7,2) - t40 * mrSges(3,1) + (-m(11) + t62) * t7 + t61 * (pkin(5) * t23 + t7) + (-m(10) * pkin(2) - mrSges(9,1)) * t24 - t65 * t19 + t64 * t17 + (t60 + t63) * t30) * g(3) + (-mrSges(1,2) + t61 * (t37 * t4 + r_base(2)) + t62 * (t15 * t37 + r_base(2)) + t63 * (pkin(15) * t37 + r_base(2)) + t58 * r_base(2) + t56 * t41 + t57 * t37) * g(2) + (-mrSges(1,1) + t62 * (t15 * t41 + r_base(1)) + t61 * (t41 * t4 + r_base(1)) + t63 * (pkin(15) * t41 + r_base(1)) + t58 * r_base(1) + t57 * t41 - t56 * t37) * g(1);
U = t1;
