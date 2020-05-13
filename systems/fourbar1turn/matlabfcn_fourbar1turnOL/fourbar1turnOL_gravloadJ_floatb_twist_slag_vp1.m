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
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:41
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnOL_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnOL_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 19:40:48
% EndTime: 2020-04-12 19:40:49
% DurationCPUTime: 0.18s
% Computational Cost: add. (67->36), mult. (113->62), div. (0->0), fcn. (85->8), ass. (0->20)
t8 = qJ(2) + qJ(3);
t6 = sin(t8);
t7 = cos(t8);
t33 = -rSges(4,1) * t7 + t6 * rSges(4,2);
t32 = rSges(4,1) * t6 + rSges(4,2) * t7;
t13 = cos(qJ(2));
t31 = pkin(2) * t13 + t33;
t11 = sin(qJ(1));
t14 = cos(qJ(1));
t30 = g(1) * t14 + g(2) * t11;
t29 = t32 * t11;
t28 = t32 * t14;
t10 = sin(qJ(2));
t22 = pkin(2) * t10;
t12 = cos(qJ(4));
t9 = sin(qJ(4));
t19 = rSges(5,1) * t12 - rSges(5,2) * t9;
t18 = rSges(3,1) * t13 - rSges(3,2) * t10;
t16 = pkin(1) + t19;
t1 = [-m(2) * (g(1) * (-t11 * rSges(2,1) - rSges(2,2) * t14) + g(2) * (rSges(2,1) * t14 - t11 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t14 - t18 * t11) + g(2) * (t11 * rSges(3,3) + t18 * t14)) - m(4) * ((g(1) * rSges(4,3) + g(2) * t31) * t14 + (g(2) * rSges(4,3) - g(1) * t31) * t11) - m(5) * ((g(1) * rSges(5,3) + g(2) * t16) * t14 + (g(2) * rSges(5,3) - g(1) * t16) * t11), -m(3) * (g(3) * t18 + t30 * (-rSges(3,1) * t10 - rSges(3,2) * t13)) - m(4) * (g(1) * (-t14 * t22 + t28) + g(2) * (-t11 * t22 + t29) + g(3) * t31), -m(4) * (g(1) * t28 + g(2) * t29 + g(3) * t33), -m(5) * (g(3) * t19 + t30 * (-rSges(5,1) * t9 - rSges(5,2) * t12)), 0];
taug = t1(:);
