% Calculate Gravitation load on the joints for
% fourbar1turnIC
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
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 11:33
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1turnIC_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fourbar1turnIC_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 11:33:32
% EndTime: 2020-05-07 11:33:33
% DurationCPUTime: 0.24s
% Computational Cost: add. (75->39), mult. (118->68), div. (4->3), fcn. (90->11), ass. (0->23)
t10 = qJ(2) + qJ(3);
t8 = sin(t10);
t9 = cos(t10);
t36 = -t9 * rSges(4,1) + t8 * rSges(4,2);
t35 = rSges(4,1) * t8 + rSges(4,2) * t9;
t15 = cos(qJ(2));
t34 = t15 * pkin(2) + t36;
t13 = sin(qJ(1));
t16 = cos(qJ(1));
t33 = g(1) * t16 + g(2) * t13;
t32 = t35 * t13;
t31 = t35 * t16;
t12 = sin(qJ(2));
t26 = pkin(2) * t12;
t24 = -qJ(4) + qJ(2);
t22 = t15 * rSges(3,1) - t12 * rSges(3,2);
t11 = sin(qJ(4));
t14 = cos(qJ(4));
t20 = t14 * rSges(5,1) - t11 * rSges(5,2);
t18 = pkin(1) + t20;
t7 = sin(qJ(3) + t24);
t5 = 0.1e1 / t7;
t1 = [-m(2) * (g(1) * (-t13 * rSges(2,1) - t16 * rSges(2,2)) + g(2) * (t16 * rSges(2,1) - t13 * rSges(2,2))) - m(3) * (g(1) * (t16 * rSges(3,3) - t22 * t13) + g(2) * (t13 * rSges(3,3) + t22 * t16)) - m(4) * ((g(1) * rSges(4,3) + g(2) * t34) * t16 + (g(2) * rSges(4,3) - g(1) * t34) * t13) - m(5) * ((g(1) * rSges(5,3) + g(2) * t18) * t16 + (g(2) * rSges(5,3) - g(1) * t18) * t13); -m(3) * (g(3) * t22 + t33 * (-rSges(3,1) * t12 - rSges(3,2) * t15)) - sin(qJ(3)) * pkin(2) / pkin(4) * t5 * m(5) * (g(3) * t20 + t33 * (-rSges(5,1) * t11 - rSges(5,2) * t14)) + (-g(1) * (-t16 * t26 + t31) - g(2) * (-t13 * t26 + t32) - g(3) * t34 - (-pkin(3) * t7 + pkin(2) * sin(t24)) / pkin(3) * t5 * (g(1) * t31 + g(2) * t32 + g(3) * t36)) * m(4);];
taug = t1(:);
