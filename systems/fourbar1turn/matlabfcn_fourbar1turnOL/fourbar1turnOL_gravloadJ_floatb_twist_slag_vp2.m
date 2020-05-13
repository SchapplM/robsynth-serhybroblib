% Calculate Gravitation load on the joints for
% fourbar1turnOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% m [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnOL_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:48
% EndTime: 2020-04-12 19:40:49
% DurationCPUTime: 0.14s
% Computational Cost: add. (65->24), mult. (112->34), div. (0->0), fcn. (85->8), ass. (0->20)
t10 = sin(qJ(2));
t13 = cos(qJ(2));
t22 = m(4) * pkin(2) + mrSges(3,1);
t8 = qJ(2) + qJ(3);
t6 = sin(t8);
t7 = cos(t8);
t31 = -mrSges(4,1) * t7 + t6 * mrSges(4,2);
t32 = mrSges(3,2) * t10 - t22 * t13 - t31;
t12 = cos(qJ(4));
t9 = sin(qJ(4));
t17 = mrSges(5,1) * t12 - mrSges(5,2) * t9;
t30 = m(5) * pkin(1) + mrSges(2,1) + t17 - t32;
t29 = mrSges(4,1) * t6 + mrSges(4,2) * t7;
t11 = sin(qJ(1));
t27 = t29 * t11;
t14 = cos(qJ(1));
t26 = t29 * t14;
t21 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t16 = mrSges(3,2) * t13 + t22 * t10;
t1 = [(t21 * g(1) - t30 * g(2)) * t14 + (t30 * g(1) + t21 * g(2)) * t11, t32 * g(3) + (t16 * t11 - t27) * g(2) + (t16 * t14 - t26) * g(1), -g(1) * t26 - g(2) * t27 - g(3) * t31, -g(3) * t17 + (g(1) * t14 + g(2) * t11) * (mrSges(5,1) * t9 + mrSges(5,2) * t12), 0];
taug = t1(:);
