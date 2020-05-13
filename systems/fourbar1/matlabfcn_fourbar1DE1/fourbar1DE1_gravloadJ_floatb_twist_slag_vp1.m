% Calculate Gravitation load on the joints for
% fourbar1DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
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
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1DE1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1DE1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1DE1_gravloadJ_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1DE1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:57:10
% EndTime: 2020-04-24 19:57:12
% DurationCPUTime: 0.29s
% Computational Cost: add. (325->53), mult. (503->86), div. (16->5), fcn. (140->4), ass. (0->33)
t45 = g(1) / 0.2e1;
t44 = g(2) / 0.2e1;
t25 = cos(qJ(1));
t40 = pkin(2) * t25;
t36 = (-0.2e1 * t40 + pkin(1)) * pkin(1);
t42 = -pkin(3) - pkin(4);
t4 = (pkin(2) - t42) * (pkin(2) + t42) + t36;
t41 = -pkin(3) + pkin(4);
t5 = (pkin(2) - t41) * (pkin(2) + t41) + t36;
t32 = sqrt(-t4 * t5);
t24 = sin(qJ(1));
t39 = t24 * pkin(2);
t34 = pkin(1) * t39;
t43 = (-t4 - t5) * t34 / t32;
t38 = rSges(3,2) * t24;
t37 = rSges(4,2) * t24;
t14 = rSges(3,1) * t39 + rSges(3,2) * t40;
t15 = rSges(4,1) * t39 + rSges(4,2) * t40;
t35 = pkin(3) ^ 2 - pkin(4) ^ 2;
t18 = pkin(2) ^ 2 + t36;
t33 = pkin(1) - t40;
t29 = 0.1e1 / pkin(3);
t17 = 0.1e1 / t18 ^ 2;
t16 = 0.1e1 / t18;
t13 = (rSges(3,1) * t25 - t38) * pkin(2);
t12 = (rSges(4,1) * t25 - t37) * pkin(2);
t11 = t18 - t35;
t10 = t18 + t35;
t9 = rSges(3,1) * t33 + pkin(2) * t38;
t8 = -rSges(3,2) * pkin(1) + t14;
t7 = -rSges(4,2) * pkin(1) + t15;
t6 = rSges(4,1) * t33 + pkin(2) * t37;
t1 = [-m(2) * (-(rSges(2,1) * g(1) + rSges(2,2) * g(2)) * t24 + (rSges(2,1) * g(2) - rSges(2,2) * g(1)) * t25) - m(3) * (((t10 * t14 + t13 * t32 + t8 * t43) * t45 + (-t10 * t13 + t14 * t32 + t9 * t43) * t44) * t29 * t16 + (g(2) * t25 + (-g(1) + (g(1) * (t9 * t16 - (t10 * t9 + t32 * t8) * t17) + g(2) * (-t8 * t16 - (-t10 * t8 + t32 * t9) * t17)) * t29 * pkin(1)) * t24) * pkin(2)) - m(4) * (((-t11 * t15 + t12 * t32 + t7 * t43) * t45 + (t11 * t12 + t15 * t32 + t6 * t43) * t44) * t16 + (g(1) * (-t6 * t16 - (-t11 * t6 + t32 * t7) * t17) + g(2) * (t7 * t16 - (t11 * t7 + t32 * t6) * t17)) * t34) / pkin(4)];
taug = t1(:);
