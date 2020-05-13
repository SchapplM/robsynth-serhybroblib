% Calculate vector of inverse dynamics joint torques for
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% qJDD [1x1]
%   Generalized joint accelerations
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
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbarprisTE_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_invdynJ_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisTE_invdynJ_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbarprisTE_invdynJ_fixb_slag_vp1: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisTE_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_invdynJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisTE_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisTE_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisTE_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:01:21
% EndTime: 2020-05-07 09:01:38
% DurationCPUTime: 7.69s
% Computational Cost: add. (5288->151), mult. (3405->330), div. (359->18), fcn. (348->2), ass. (0->127)
t64 = (qJ(1) + pkin(3));
t137 = -pkin(2) + t64;
t49 = pkin(1) + t137;
t138 = -pkin(2) - t64;
t50 = pkin(1) + t138;
t159 = t49 * t50;
t48 = pkin(1) - t137;
t134 = t48 * t159;
t47 = pkin(1) - t138;
t116 = t47 * t134;
t76 = sqrt(-t116);
t20 = 0.1e1 / t76;
t59 = 0.1e1 / t64 ^ 2;
t157 = t59 * t20;
t10 = -t134 + (t159 + (t49 - t50) * t48) * t47;
t68 = (qJ(1) ^ 2);
t70 = (pkin(3) ^ 2);
t185 = (t70 + t68);
t110 = -2 * pkin(3) * qJ(1) - t185;
t72 = (pkin(2) ^ 2);
t73 = (pkin(1) ^ 2);
t149 = (t73 - t72);
t27 = t110 + t149;
t58 = 1 / t64;
t163 = t27 * t58;
t21 = t20 / t116;
t114 = t21 * t163 / 0.2e1;
t162 = t27 * t59;
t175 = -2 * t64;
t168 = (t10 * t114 + ((t58 * t175) - t162) * t20) * qJD(1) * t58;
t189 = -t168 / 0.2e1;
t142 = t27 * qJD(1);
t131 = t20 * t142;
t113 = -t131 / 0.2e1;
t104 = t59 * t113;
t57 = t64 ^ 2;
t75 = 1 / (pkin(1) ^ 2);
t187 = 1 / t57 ^ 2 * t75;
t35 = (m(2) * rSges(2,2) - m(3) * rSges(3,1));
t36 = (m(2) * rSges(2,1) + m(3) * rSges(3,3));
t24 = (-t35 * g(1) + g(2) * t36);
t184 = t64 * t157;
t179 = 0.1e1 / t50;
t180 = 0.1e1 / t48;
t183 = t179 * t180;
t37 = 0.1e1 / t47;
t41 = 0.1e1 / t49;
t160 = t37 * t41;
t117 = t160 * t183;
t74 = 0.1e1 / pkin(1);
t178 = t74 ^ 2;
t40 = 0.1e1 / t48 ^ 2;
t44 = 0.1e1 / t50 ^ 2;
t33 = t64 * qJD(1);
t177 = -0.2e1 * t33;
t176 = 0.2e1 * t33;
t67 = qJ(1) * t68;
t174 = 2 * t67;
t173 = 5 * t68;
t172 = (m(3) * g(2));
t171 = -t20 / 0.2e1;
t170 = -t58 / 0.2e1;
t169 = t59 / 0.2e1;
t9 = qJD(1) * t10;
t167 = t21 * t9;
t166 = qJ(1) / 0.2e1;
t26 = -t110 + t149;
t14 = (rSges(3,1) * t26) - rSges(3,3) * t76;
t165 = t20 * t14;
t164 = (t24 * t73);
t161 = t33 * t59;
t158 = (pkin(1) + t64) * (pkin(1) - t64);
t83 = t64 * t57;
t60 = 1 / t83;
t156 = t60 * t20;
t66 = qJD(1) ^ 2;
t155 = t64 * t66;
t154 = t66 * t59;
t153 = t66 * t60;
t148 = qJ(1) * t70;
t147 = qJD(1) * t9;
t146 = t26 * qJ(1);
t145 = -qJDD(1) / 0.2e1;
t144 = qJD(1) * t33;
t143 = qJD(1) * t74;
t141 = t60 * qJD(1);
t140 = 2 * m(4) * t83;
t139 = 0.2e1 * t161;
t136 = 2 * qJ(1);
t135 = t27 * t157;
t127 = rSges(3,3) * t171;
t133 = (rSges(3,1) * t176 + t127 * t9) * t171;
t132 = -t58 * t9 / 0.4e1;
t130 = rSges(2,2) * t20 / 0.2e1;
t129 = t10 * t171;
t128 = rSges(3,1) * t171;
t125 = -t156 / 0.2e1;
t124 = t154 / 0.2e1;
t123 = t146 / 0.2e1;
t120 = t142 / 0.2e1;
t107 = t166 - t165 / 0.2e1;
t16 = -rSges(3,1) * t76 - rSges(3,3) * t26;
t106 = t123 - t16 / 0.2e1;
t105 = qJD(1) * t117;
t103 = t120 * t157;
t25 = t27 ^ 2;
t100 = t25 * t66 * t117;
t99 = -qJ(1) * t117 / 0.4e1;
t38 = 0.1e1 / t47 ^ 2;
t42 = 0.1e1 / t49 ^ 2;
t97 = t25 * t38 * t40 * t42 * t44 * t147;
t96 = t105 * t187;
t1 = t114 * t147 + (-0.2e1 * t58 * t144 + (qJDD(1) * t58 - t154) * t27) * t20;
t94 = t1 * t58 / 0.2e1 + t27 * t66 * t125;
t88 = t107 * t162 + t170 * t26;
t87 = -t106 * t135 + t170 * t76;
t77 = (pkin(2) * t72);
t69 = pkin(3) * t70;
t63 = t69 * t172;
t17 = -rSges(2,1) * t76 - rSges(2,2) * t26;
t15 = -rSges(2,1) * t26 + rSges(2,2) * t76;
t8 = rSges(2,1) * t171 * t9 + rSges(2,2) * t177;
t7 = rSges(3,3) * t177 + t128 * t9;
t6 = rSges(2,1) * t177 + t130 * t9;
t4 = t88 * t143;
t3 = t87 * t143;
t2 = [m(2) * (-t17 * ((rSges(2,1) * t129 + rSges(2,2) * t175) * t104 + (t120 * t156 + t189) * t17) + ((rSges(2,1) * t175 + t10 * t130) * t103 + (t168 / 0.2e1 + t60 * t113) * t15) * t15) * t178 * t104 + m(2) * (-(t15 * t6 + t17 * t8) * t25 * t96 / 0.2e1 + (-(t104 * t8 - t17 * t94) * t17 + (t103 * t6 + t15 * t94) * t15) * t178 * t135 + (t15 ^ 2 + t17 ^ 2) * (t27 * t33 * t96 + (t58 * t100 - t97 / 0.4e1) * t187)) / 0.2e1 + m(3) * ((t26 * t124 + (t26 * t145 - t144 + (t76 * t166 - t14 / 0.2e1) * t1) * t58 + (-t107 * t153 + (t66 / 0.2e1 + (t9 * t99 + t133) * qJD(1)) * t59) * t27) * t88 + (t76 * t124 + (-t1 * t106 + t145 * t76) * t58 + (qJD(1) * t132 + (t106 * t153 + (-t66 * t26 / 0.2e1 + (-t33 * qJ(1) + t7 / 0.2e1) * qJD(1)) * t59) * t27) * t20) * t87) * t178 + (Icges(3,2) + Icges(2,3)) * (t1 * t20 * t163 - t131 * t168 + t27 * t105 * t139 - t59 * t97 / 0.2e1 + t60 * t100) + (((((rSges(4,1) * g(1) + rSges(4,2) * g(2)) * t140) + (-pkin(3) * t77 + (t69 + 4 * t148 + (t173 + t73) * pkin(3) + t174) * pkin(2)) * g(1) * m(3) + ((t77 + (pkin(3) * t136 + t185 - t73) * pkin(2)) * (t36 * g(1) + t35 * g(2)))) * t76 + ((-pkin(3) * t172 + t24) * t72 * t77) + (2 * (-t164 + t63 + (3 * t148 + (3 * t68 + t73) * pkin(3) + t67) * t172) * t77) + ((t63 + (4 * qJ(1) * t172 + t24) * t70 + (t24 * t136 + (t173 - t73) * t172) * pkin(3) + t172 * t174 + t24 * t68 + t164) * pkin(2) * t158) + ((t72 + t158) * (rSges(4,1) * g(2) - rSges(4,2) * g(1)) * t140)) * t20 / pkin(2) * t169 + (-((t129 * t3 + t175 * t4) * t58 + 0.2e1 * (t4 * (-(t60 * qJ(1)) / 0.2e1 - rSges(3,1) * t184 + t14 * t156 / 0.2e1 + (t59 * t99 - t127 * t157 / 0.2e1) * t10) + t3 * (-qJ(1) * t184 + t123 * t156 + (rSges(3,3) * t175 + t10 * t128) * t157 / 0.2e1 + t16 * t125)) * t27) * qJD(1) / 0.2e1 + (-t14 * t4 + t16 * t3 + (-t26 * t3 + t4 * t76) * qJ(1)) * t189 + t4 * ((-t107 * t59 + t170) * t176 + ((-t14 * t167 / 0.4e1 + t133) * t59 + (-qJ(1) + t165) * t141) * t27) + t3 * (-(t146 / 0.4e1 - t16 / 0.4e1) * t162 * t167 + (t132 + t106 * t139 + (-t16 * t141 + t7 * t169 + (t141 * t26 - t161) * qJ(1)) * t27) * t20)) * m(3)) * t74 + ((0.4e1 * t155 + (-0.32e2 * qJDD(1) * t57 - 0.64e2 * t155) * t73 * t72 * t75 / (pkin(2) ^ 2) / 0.8e1) * t117 + (-0.2e1 * (t179 * t40 + t180 * t44) * t160 + 0.2e1 * (t37 * t42 + t38 * t41) * t183) * t57 * t66) * ((rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3));];
tau = t2;
