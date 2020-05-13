% Calculate Gravitation load on the joints for
% palh2m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m2IC_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2IC_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2IC_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2IC_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2IC_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2IC_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 06:51:27
% EndTime: 2020-05-03 06:51:27
% DurationCPUTime: 0.53s
% Computational Cost: add. (402->60), mult. (402->82), div. (0->0), fcn. (266->13), ass. (0->38)
t83 = m(6) + m(7);
t32 = qJ(2) + qJ(3);
t29 = qJ(4) + t32;
t24 = cos(t29);
t39 = cos(qJ(2));
t31 = t39 * pkin(4);
t7 = pkin(2) * cos(t32) + t31 + pkin(1);
t82 = pkin(5) * t24 + t7;
t35 = sin(qJ(2));
t38 = cos(qJ(3));
t26 = rSges(4,1) * m(4) * t38;
t34 = sin(qJ(3));
t67 = rSges(4,2) * t34;
t4 = m(3) * rSges(3,1) + t26 + (pkin(4) - t67) * m(4);
t9 = (rSges(4,1) * t34 + rSges(4,2) * t38) * m(4);
t6 = rSges(3,2) * m(3) + t9;
t81 = t6 * t35 - t4 * t39;
t36 = sin(qJ(1));
t40 = cos(qJ(1));
t80 = -g(1) * t40 - t36 * g(2);
t79 = m(5) + t83;
t72 = m(5) * rSges(5,2);
t61 = -m(2) * rSges(2,2) + m(3) * rSges(3,3) + rSges(4,3) * m(4);
t33 = sin(qJ(6));
t37 = cos(qJ(6));
t60 = t37 * rSges(7,1) - rSges(7,2) * t33;
t12 = rSges(7,1) * t33 + rSges(7,2) * t37;
t23 = sin(t29);
t50 = t24 * rSges(5,1) - rSges(5,2) * t23 + t7;
t25 = qJ(5) + t29;
t17 = sin(t25);
t18 = cos(t25);
t48 = t18 * rSges(6,1) - rSges(6,2) * t17 + t82;
t44 = -m(2) * rSges(2,1) + (-m(3) - m(4)) * pkin(1) + t81;
t42 = -rSges(7,3) * t17 + (pkin(3) + t60) * t18 + t82;
t14 = t83 * pkin(5) + m(5) * rSges(5,1);
t10 = -m(4) * t67 + t26;
t1 = [-(t44 * g(1) + t61 * g(2)) * t36 + (-t61 * g(1) + t44 * g(2)) * t40 - m(5) * ((g(1) * rSges(5,3) + g(2) * t50) * t40 + (g(2) * rSges(5,3) - g(1) * t50) * t36) - m(6) * ((g(1) * rSges(6,3) + g(2) * t48) * t40 + (g(2) * rSges(6,3) - g(1) * t48) * t36) - m(7) * ((-g(1) * t42 - g(2) * t12) * t36 + (-g(1) * t12 + g(2) * t42) * t40); (t10 * t39 - t79 * t31 - t9 * t35 + t81) * g(3) + ((-t6 + t9) * t39 + (-t79 * pkin(4) + t10 - t4) * t35) * t80; (g(3) * t72 - t14 * t80) * t23 - (g(3) * t14 + t72 * t80) * t24; (-t60 * (-t36 * g(1) + g(2) * t40) + (g(3) * t17 - t18 * t80) * t12) * m(7);];
taug = t1(:);
