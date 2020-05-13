% Calculate Gravitation load on the joints for
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% m [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbarprisTE_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_gravloadJ_floatb_twist_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisTE_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_gravloadJ_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisTE_gravloadJ_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisTE_gravloadJ_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:01:21
% EndTime: 2020-05-07 09:01:23
% DurationCPUTime: 0.26s
% Computational Cost: add. (104->46), mult. (207->64), div. (3->4), fcn. (4->2), ass. (0->25)
t17 = (qJ(1) + pkin(3));
t47 = 2 * t17 ^ 2;
t7 = (m(2) * rSges(2,2) - m(3) * rSges(3,1));
t8 = (m(2) * rSges(2,1) + m(3) * rSges(3,3));
t1 = (-t7 * g(1) + g(2) * t8);
t20 = (qJ(1) ^ 2);
t19 = (qJ(1) * t20);
t45 = 2 * t19;
t44 = 5 * t20;
t43 = (m(3) * g(2));
t24 = (pkin(1) ^ 2);
t41 = (t1 * t24);
t22 = (pkin(3) ^ 2);
t39 = (qJ(1) * t22);
t38 = m(4) * t17 * t47;
t37 = -pkin(2) - t17;
t36 = -pkin(2) + t17;
t35 = 2 * qJ(1);
t33 = (rSges(4,1) * g(2) - rSges(4,2) * g(1)) * t38;
t26 = pkin(2) ^ 2;
t27 = pkin(2) * t26;
t25 = sqrt(-((pkin(1) + t37) * (pkin(1) + t36) * (pkin(1) - t36) * (pkin(1) - t37)));
t21 = (pkin(3) * t22);
t16 = (t21 * t43);
t2 = [(((rSges(4,1) * g(1) + rSges(4,2) * g(2)) * t38 + (-pkin(3) * t27 + (t21 + 4 * t39 + (t44 + t24) * pkin(3) + t45) * pkin(2)) * g(1) * m(3) + (t27 + (pkin(3) * t35 + t20 + t22 - t24) * pkin(2)) * (t8 * g(1) + t7 * g(2))) * t25 + (2 * (-t41 + t16 + (3 * t39 + (3 * t20 + t24) * pkin(3) + t19) * t43) * t27) + (((t16 + (4 * qJ(1) * t43 + t1) * t22 + (t1 * t35 + (t44 - t24) * t43) * pkin(3) + t43 * t45 + t1 * t20 + t41) * pkin(2) + t33) * (pkin(1) - t17) * (pkin(1) + t17)) + (((-pkin(3) * t43 + t1) * t27 + t33) * t26)) / t25 / pkin(1) / pkin(2) / t47];
taug = t2(:);
