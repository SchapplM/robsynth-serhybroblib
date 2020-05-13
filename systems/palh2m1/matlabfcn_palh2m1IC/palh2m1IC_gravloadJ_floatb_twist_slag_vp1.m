% Calculate Gravitation load on the joints for
% palh2m1IC
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
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:04
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m1IC_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1IC_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1IC_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1IC_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1IC_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1IC_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:03:48
% EndTime: 2020-05-03 01:03:48
% DurationCPUTime: 0.40s
% Computational Cost: add. (325->79), mult. (393->106), div. (0->0), fcn. (238->12), ass. (0->49)
t31 = sin(qJ(5));
t35 = cos(qJ(5));
t19 = t35 * rSges(6,1) - t31 * rSges(6,2);
t61 = -pkin(4) - t19;
t30 = qJ(2) + qJ(3);
t29 = qJ(4) + t30;
t24 = cos(t29);
t77 = t61 * t24;
t73 = rSges(6,3) + pkin(6);
t75 = t24 * t73;
t38 = cos(qJ(1));
t74 = t38 * g(2);
t32 = sin(qJ(3));
t36 = cos(qJ(3));
t10 = (rSges(4,1) * t32 + rSges(4,2) * t36) * m(4);
t23 = sin(t29);
t47 = t24 * rSges(5,1) - t23 * rSges(5,2);
t33 = sin(qJ(2));
t37 = cos(qJ(2));
t26 = rSges(4,1) * m(4) * t36;
t60 = rSges(4,2) * t32;
t4 = rSges(3,1) * m(3) + t26 + (pkin(2) - t60) * m(4);
t8 = rSges(3,2) * m(3) + t10;
t72 = -t8 * t33 + t4 * t37;
t34 = sin(qJ(1));
t51 = t38 * g(1) + t34 * g(2);
t66 = t37 * pkin(2);
t69 = pkin(3) * cos(t30);
t9 = pkin(1) + t66 + t69;
t71 = t9 * t74;
t70 = pkin(3) * sin(t30);
t54 = t34 * t75;
t53 = t38 * t75;
t52 = m(2) * rSges(2,2) + m(3) * rSges(3,3) + rSges(4,3) * m(4);
t20 = rSges(5,2) * m(5) - m(6) * t73;
t3 = m(5) * rSges(5,1) - m(6) * t61;
t50 = -(-g(3) * t20 + t51 * t3) * t23 - (g(3) * t3 + t51 * t20) * t24;
t46 = -rSges(5,1) * t23 - rSges(5,2) * t24;
t44 = -t69 - t47;
t42 = -m(2) * rSges(2,1) - t72 + (-m(3) - m(4)) * pkin(1);
t41 = -t69 + t77;
t40 = t73 * t23 - t77;
t39 = (-g(3) * t73 + t51 * t61) * t23;
t18 = t31 * rSges(6,1) + t35 * rSges(6,2);
t12 = -t33 * pkin(2) - t70;
t11 = -m(4) * t60 + t26;
t7 = t38 * t12;
t6 = t34 * t12;
t1 = [-(t42 * g(1) - t52 * g(2)) * t34 + (t52 * g(1) + t42 * g(2)) * t38 - m(5) * (t71 + (-g(1) * rSges(5,3) + g(2) * t47) * t38 + (g(1) * (-t47 - t9) - g(2) * rSges(5,3)) * t34) - m(6) * (t71 + (-g(1) * t18 + g(2) * t40) * t38 + (g(1) * (-t40 - t9) - g(2) * t18) * t34); -m(5) * (g(1) * (t46 * t38 + t7) + g(2) * (t46 * t34 + t6)) - m(6) * (g(1) * (t7 + t53) + g(2) * (t6 + t54) + t39) + t50 + t51 * (t33 * t4 + t37 * t8) + (-m(5) * (t44 - t66) - m(6) * (t41 - t66) + t72) * g(3); -m(6) * (g(1) * (-t38 * t70 + t53) + g(2) * (-t34 * t70 + t54) + t39) + t50 + (-m(5) * t44 - m(6) * t41 - t10 * t33 + t11 * t37) * g(3) + (t10 * t37 + t11 * t33 - m(5) * (t46 - t70)) * t51; (-t19 * (-t34 * g(1) + t74) + (-g(3) * t23 + t51 * t24) * t18) * m(6);];
taug = t1(:);
