% Calculate Gravitation load on the joints for
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m2DE2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2DE2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [22x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2DE2_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2DE2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:58:29
% EndTime: 2020-05-02 20:58:32
% DurationCPUTime: 0.54s
% Computational Cost: add. (489->93), mult. (542->116), div. (0->0), fcn. (404->35), ass. (0->56)
t51 = sin(qJ(3));
t57 = cos(qJ(3));
t88 = m(5) + m(6);
t69 = pkin(5) * t88 + m(4) * rSges(4,1);
t87 = m(4) * rSges(4,2);
t23 = -t51 * t87 + t69 * t57;
t22 = t69 * t51 + t57 * t87;
t55 = sin(pkin(17));
t61 = cos(pkin(17));
t54 = sin(pkin(18));
t60 = cos(pkin(18));
t72 = t60 * rSges(7,1) + rSges(7,2) * t54;
t76 = m(11) + m(4) + m(8) + t88;
t93 = rSges(7,1) * t54 - t60 * rSges(7,2);
t6 = t76 * pkin(1) + rSges(3,1) * m(3) + (t55 * t93 + t61 * t72) * m(7) + t22;
t92 = rSges(3,2) * m(3) + (t72 * t55 - t93 * t61) * m(7) - t23;
t50 = sin(qJ(4));
t56 = cos(qJ(4));
t73 = t50 * rSges(6,1) + t56 * rSges(6,2);
t91 = -m(2) * rSges(2,2) + m(11) * rSges(11,3) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3) + t73 * m(6) + m(7) * rSges(7,3) + m(8) * rSges(8,3) + m(9) * rSges(9,3) + m(10) * rSges(10,3);
t52 = sin(qJ(2));
t58 = cos(qJ(2));
t90 = -(m(10) + m(9) + m(3) + t76) * pkin(15) - m(2) * rSges(2,1) + pkin(14) * m(7) + t6 * t52 + t92 * t58;
t85 = m(9) * rSges(9,2);
t48 = sin(pkin(19));
t49 = cos(pkin(19));
t27 = t57 * t48 + t51 * t49;
t28 = t51 * t48 - t57 * t49;
t18 = qJ(2) + atan2(t28, t27);
t47 = pkin(18) - pkin(22);
t46 = -qJ(2) + t47;
t80 = pkin(21) - atan2(cos(t46), -sin(t46));
t78 = -qJ(2) + t80;
t77 = -pkin(21) + t47;
t53 = sin(qJ(1));
t59 = cos(qJ(1));
t35 = g(1) * t59 + g(2) * t53;
t75 = g(1) * t53 - g(2) * t59;
t74 = t56 * rSges(6,1) - t50 * rSges(6,2);
t40 = -qJ(3) - qJ(2) + t77;
t21 = atan2(-sin(t40), cos(t40));
t14 = -t21 + t78;
t11 = sin(t14);
t12 = cos(t14);
t15 = -t21 + t80;
t16 = sin(t18);
t17 = cos(t18);
t44 = pkin(2) * m(10) + m(9) * rSges(9,1);
t64 = (g(3) * t85 + t35 * t44) * t17 + (g(3) * t44 - t35 * t85) * t16 + ((-rSges(11,1) * t12 - rSges(11,2) * t11) * t35 + ((t58 * rSges(11,1) - t52 * rSges(11,2)) * sin(t15) - (t52 * rSges(11,1) + rSges(11,2) * t58) * cos(t15)) * g(3)) * m(11);
t43 = -pkin(20) + t77;
t38 = cos(t43);
t37 = sin(t43);
t4 = atan2(t28, -t27) + t18;
t3 = cos(t4);
t2 = sin(t4);
t1 = [(-t91 * g(1) + t90 * g(2)) * t59 - (t90 * g(1) + t91 * g(2)) * t53 + ((-cos(t47) * rSges(8,1) + sin(t47) * rSges(8,2)) * m(8) + (rSges(10,1) * t2 + rSges(10,2) * t3) * m(10) + (pkin(4) * sin(t78) - rSges(11,1) * t11 + rSges(11,2) * t12) * m(11) - t16 * t44 - t17 * t85 - t38 * (m(5) * rSges(5,1) + (pkin(9) + t74) * m(6)) + t37 * ((-rSges(6,3) - pkin(11)) * m(6) + m(5) * rSges(5,2))) * t75, (g(3) * t92 + t35 * t6) * t58 + (g(3) * t6 - t35 * t92) * t52 + (-(t35 * rSges(10,1) + g(3) * rSges(10,2)) * t3 + (-g(3) * rSges(10,1) + t35 * rSges(10,2)) * t2) * m(10) + t64, (-t23 * g(3) + t35 * t22) * t58 + (t22 * g(3) + t35 * t23) * t52 + t64, -(t74 * t75 + (g(3) * t37 + t35 * t38) * t73) * m(6)];
taug = t1(:);
