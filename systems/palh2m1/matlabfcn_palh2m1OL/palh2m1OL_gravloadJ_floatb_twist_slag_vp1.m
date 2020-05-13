% Calculate Gravitation load on the joints for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m1OL_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1OL_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:07:25
% EndTime: 2020-05-03 00:07:28
% DurationCPUTime: 0.48s
% Computational Cost: add. (299->77), mult. (357->106), div. (0->0), fcn. (226->12), ass. (0->48)
t29 = sin(qJ(5));
t33 = cos(qJ(5));
t73 = -rSges(6,1) * t33 + rSges(6,2) * t29;
t58 = -pkin(4) + t73;
t28 = qJ(2) + qJ(3);
t27 = qJ(4) + t28;
t22 = cos(t27);
t74 = t58 * t22;
t70 = rSges(6,3) + pkin(6);
t72 = t22 * t70;
t36 = cos(qJ(1));
t71 = g(2) * t36;
t30 = sin(qJ(3));
t34 = cos(qJ(3));
t8 = (rSges(4,1) * t30 + rSges(4,2) * t34) * m(4);
t21 = sin(t27);
t44 = rSges(5,1) * t22 - t21 * rSges(5,2);
t32 = sin(qJ(1));
t48 = g(1) * t36 + t32 * g(2);
t24 = rSges(4,1) * m(4) * t34;
t55 = rSges(4,2) * t30;
t2 = rSges(3,1) * m(3) + t24 + (pkin(2) - t55) * m(4);
t31 = sin(qJ(2));
t35 = cos(qJ(2));
t6 = rSges(3,2) * m(3) + t8;
t69 = t2 * t35 - t6 * t31;
t65 = pkin(3) * cos(t28);
t67 = pkin(2) * t35;
t7 = pkin(1) + t65 + t67;
t68 = t7 * t71;
t66 = pkin(3) * sin(t28);
t51 = t32 * t72;
t50 = t36 * t72;
t49 = m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3);
t43 = -rSges(5,1) * t21 - rSges(5,2) * t22;
t42 = -t65 - t44;
t40 = -m(2) * rSges(2,1) - t69 + (-m(3) - m(4)) * pkin(1);
t39 = -t65 + t74;
t38 = t70 * t21 - t74;
t37 = (-g(3) * t70 + t48 * t58) * t21;
t18 = rSges(5,2) * m(5) - m(6) * t70;
t16 = rSges(6,1) * t29 + rSges(6,2) * t33;
t10 = -pkin(2) * t31 - t66;
t9 = -m(4) * t55 + t24;
t5 = t36 * t10;
t4 = t32 * t10;
t1 = m(5) * rSges(5,1) - m(6) * t58;
t3 = [-(t40 * g(1) - t49 * g(2)) * t32 + (t49 * g(1) + t40 * g(2)) * t36 - m(5) * (t68 + (-g(1) * rSges(5,3) + g(2) * t44) * t36 + (g(1) * (-t44 - t7) - g(2) * rSges(5,3)) * t32) - m(6) * (t68 + (-g(1) * t16 + g(2) * t38) * t36 + (g(1) * (-t38 - t7) - g(2) * t16) * t32), -m(5) * (g(1) * (t43 * t36 + t5) + g(2) * (t43 * t32 + t4)) - m(6) * (g(1) * (t5 + t50) + g(2) * (t4 + t51) + t37) + t48 * (t2 * t31 + t35 * t6) + (-m(5) * (t42 - t67) - m(6) * (t39 - t67) + t69) * g(3), -m(6) * (g(1) * (-t36 * t66 + t50) + g(2) * (-t32 * t66 + t51) + t37) + (-m(5) * t42 - m(6) * t39 - t8 * t31 + t9 * t35) * g(3) + (t31 * t9 + t35 * t8 - m(5) * (t43 - t66)) * t48, (-g(3) * t18 + t48 * t1) * t21 + (g(3) * t1 + t48 * t18) * t22, (t73 * (-t32 * g(1) + t71) + (-g(3) * t21 + t48 * t22) * t16) * m(6)];
taug = t3(:);
