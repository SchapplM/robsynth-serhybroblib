% Calculate potential energy for
% palh1m2TE
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
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m2TE_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m2TE_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2TE_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_energypot_floatb_twist_slag_vp2: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2TE_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2TE_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 19:58:53
% EndTime: 2020-05-01 19:58:54
% DurationCPUTime: 0.33s
% Computational Cost: add. (404->79), mult. (418->77), div. (0->0), fcn. (308->20), ass. (0->46)
t55 = m(5) + m(6);
t46 = m(11) + m(4) + m(8) + t55;
t54 = -m(3) - m(9) - m(10) - t46;
t30 = sin(qJ(4));
t35 = cos(qJ(4));
t53 = -mrSges(6,1) * t35 + mrSges(6,2) * t30;
t18 = pkin(11) * m(6) - mrSges(5,2) + mrSges(6,3);
t20 = pkin(9) * m(6) + mrSges(5,1);
t24 = sin(pkin(20));
t28 = cos(pkin(20));
t44 = t18 * t28 + t24 * t20;
t43 = mrSges(6,1) * t30 + mrSges(6,2) * t35 - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t21 = pkin(2) * m(10) + mrSges(9,1);
t25 = sin(pkin(19));
t29 = cos(pkin(19));
t42 = mrSges(9,2) * t29 + t21 * t25 + mrSges(11,2) + mrSges(4,2);
t41 = m(7) + m(2) - t54;
t22 = sin(pkin(22));
t26 = cos(pkin(22));
t11 = m(11) * pkin(4) - t18 * t24 + t20 * t28;
t23 = sin(pkin(21));
t27 = cos(pkin(21));
t6 = t11 * t27 - t44 * t23 + mrSges(8,1);
t7 = t11 * t23 + t44 * t27 - mrSges(8,2);
t1 = -t7 * t22 + t6 * t26;
t2 = t6 * t22 + t7 * t26;
t10 = t55 * pkin(5) - t25 * mrSges(9,2) + t21 * t29 + mrSges(11,1) + mrSges(4,1);
t34 = sin(pkin(17));
t39 = cos(pkin(17));
t16 = mrSges(7,1) * t34 + mrSges(7,2) * t39;
t17 = mrSges(7,1) * t39 - mrSges(7,2) * t34;
t31 = sin(qJ(3));
t33 = sin(pkin(18));
t36 = cos(qJ(3));
t38 = cos(pkin(18));
t3 = t46 * pkin(1) + t10 * t31 + t16 * t33 + t17 * t38 + t42 * t36 + mrSges(3,1) + mrSges(10,1);
t32 = sin(qJ(2));
t37 = cos(qJ(2));
t5 = t10 * t36 - t16 * t38 + t17 * t33 - t42 * t31 - mrSges(3,2) - mrSges(10,2);
t12 = t28 * t23 + t24 * t27;
t13 = -t24 * t23 + t28 * t27;
t8 = t12 * t26 + t22 * t13;
t9 = -t22 * t12 + t13 * t26;
t40 = pkin(14) * m(7) + t1 * t38 + t2 * t33 + t3 * t32 - t5 * t37 - mrSges(2,1) + t54 * pkin(15) - t53 * (t33 * t8 + t38 * t9);
t14 = m(1) + t41;
t4 = (t40 * g(1) + t43 * g(2)) * cos(qJ(1)) + (-t43 * g(1) + t40 * g(2)) * sin(qJ(1)) - mrSges(1,1) * g(1) - mrSges(1,2) * g(2) + (-g(1) * r_base(1) - g(2) * r_base(2)) * t14 + (t53 * (-t9 * t33 + t8 * t38) + m(7) * pkin(16) - t41 * pkin(13) + t1 * t33 - t14 * r_base(3) - t2 * t38 - t3 * t37 - t5 * t32 - mrSges(1,3) - mrSges(2,3)) * g(3);
U = t4;
