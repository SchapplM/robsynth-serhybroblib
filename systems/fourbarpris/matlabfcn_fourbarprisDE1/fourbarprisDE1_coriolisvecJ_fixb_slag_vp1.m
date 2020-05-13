% Calculate vector of centrifugal and Coriolis load on the joints for
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
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
% tauc [1x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = fourbarprisDE1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_coriolisvecJ_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE1_coriolisvecJ_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_coriolisvecJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE1_coriolisvecJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisDE1_coriolisvecJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisDE1_coriolisvecJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:09:52
% EndTime: 2020-05-07 09:10:03
% DurationCPUTime: 3.27s
% Computational Cost: add. (4990->108), mult. (3063->252), div. (338->15), fcn. (326->2), ass. (0->102)
t52 = (qJ(1) + pkin(3));
t112 = -pkin(2) + t52;
t41 = pkin(1) + t112;
t113 = -pkin(2) - t52;
t42 = pkin(1) + t113;
t124 = t41 * t42;
t40 = pkin(1) - t112;
t109 = t40 * t124;
t39 = pkin(1) - t113;
t94 = t39 * t109;
t60 = sqrt(-t94);
t20 = 0.1e1 / t60;
t137 = -t20 / 0.2e1;
t10 = -t109 + (t124 + (t41 - t42) * t40) * t39;
t9 = qJD(1) * t10;
t107 = t9 * t137;
t28 = (qJD(1) * t52);
t140 = -2 * t28;
t152 = rSges(3,1) * t107 / 0.2e1 + (rSges(3,3) * t140) / 0.2e1;
t48 = 0.1e1 / t52 ^ 2;
t123 = t48 * t20;
t21 = t20 / t94;
t149 = (pkin(1) ^ 2);
t119 = (-pkin(2) ^ 2 + t149);
t89 = -qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t26 = t89 + t119;
t101 = t21 * t26 / 0.2e1;
t127 = t26 * t48;
t138 = -2 * t52;
t47 = 1 / t52;
t147 = qJD(1) * t47;
t132 = (t47 * t10 * t101 + ((t47 * t138) - t127) * t20) * t147;
t151 = -t132 / 0.2e1;
t116 = t26 * qJD(1);
t105 = t20 * t116;
t91 = -t105 / 0.2e1;
t84 = t48 * t91;
t46 = t52 ^ 2;
t148 = 1 / t46 ^ 2 / t149;
t142 = 0.1e1 / t42;
t143 = 0.1e1 / t40;
t146 = t142 * t143;
t27 = (rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3);
t144 = t27 * t146;
t58 = 1 / pkin(1);
t141 = t58 ^ 2;
t32 = 0.1e1 / t40 ^ 2;
t36 = 0.1e1 / t42 ^ 2;
t139 = 2 * t28;
t136 = -t47 / 0.2e1;
t53 = qJD(1) ^ 2;
t134 = t53 / 0.2e1;
t111 = t26 * t123;
t1 = -t53 * t111 + (t101 * t9 + t140 * t20) * t147;
t133 = t1 * t47;
t131 = t21 * t9;
t130 = qJ(1) / 0.2e1;
t25 = -t89 + t119;
t14 = rSges(3,1) * t25 - rSges(3,3) * t60;
t129 = t20 * t14;
t49 = t47 / t46;
t128 = t20 * t49;
t126 = t28 * t48;
t29 = 0.1e1 / t39;
t33 = 0.1e1 / t41;
t125 = t29 * t33;
t122 = t49 * t53;
t118 = t25 * qJ(1);
t117 = qJD(1) * t58;
t115 = t49 * qJD(1);
t114 = 0.2e1 * t126;
t108 = (rSges(3,1) * t139 + rSges(3,3) * t107) * t137;
t106 = -t47 * t9 / 0.4e1;
t104 = t20 * rSges(2,2) / 0.2e1;
t103 = t10 * t137;
t102 = -t128 / 0.2e1;
t100 = t48 * t134;
t99 = t118 / 0.2e1;
t98 = t116 / 0.2e1;
t95 = (t125 * t146);
t87 = t130 - t129 / 0.2e1;
t16 = -rSges(3,1) * t60 - (rSges(3,3) * t25);
t86 = t99 - t16 / 0.2e1;
t85 = qJD(1) * t95;
t83 = t98 * t123;
t24 = t26 ^ 2;
t80 = t24 * t53 * t95;
t79 = -(qJ(1) * t95) / 0.4e1;
t30 = 0.1e1 / t39 ^ 2;
t34 = 0.1e1 / t41 ^ 2;
t77 = qJD(1) * t24 * t30 * t32 * t34 * t36 * t9;
t76 = t85 * t148;
t74 = t133 / 0.2e1 + t26 * t53 * t102;
t68 = t127 * t87 + t136 * t25;
t67 = -t111 * t86 + t136 * t60;
t17 = -rSges(2,1) * t60 - rSges(2,2) * t25;
t15 = -rSges(2,1) * t25 + rSges(2,2) * t60;
t8 = rSges(2,1) * t107 + rSges(2,2) * t140;
t6 = rSges(2,1) * t140 + t104 * t9;
t4 = t68 * t117;
t3 = t67 * t117;
t2 = [m(2) * ((-(-t17 * t74 + t8 * t84) * t17 + (t15 * t74 + t6 * t83) * t15) * t141 * t111 - (t6 * t15 + t8 * t17) * t24 * t76 / 0.2e1 + (t15 ^ 2 + t17 ^ 2) * ((t26 * t28 * t76) + ((t47 * t80) - t77 / 0.4e1) * t148)) / 0.2e1 + m(2) * (-t17 * ((rSges(2,1) * t103 + rSges(2,2) * t138) * t84 + (t98 * t128 + t151) * t17) + ((rSges(2,1) * t138 + t10 * t104) * t83 + (t132 / 0.2e1 + t49 * t91) * t15) * t15) * t141 * t84 + 0.2e1 * (t34 * t29 + t33 * t30) * t46 * t53 * t144 + (((t25 * t100 + (-(t28 * qJD(1)) + (t60 * t130 - t14 / 0.2e1) * t1) * t47 + (-t87 * t122 + (t134 + (t79 * t9 + t108) * qJD(1)) * t48) * t26) * t68 + (t60 * t100 - t86 * t133 + (qJD(1) * t106 + (t86 * t122 + (-(t53 * t25) / 0.2e1 + (-(t28 * qJ(1)) + t152) * qJD(1)) * t48) * t26) * t20) * t67) * t141 + (t4 * ((-t48 * t87 + t136) * t139 + ((-t14 * t131 / 0.4e1 + t108) * t48 + (-qJ(1) + t129) * t115) * t26) + t3 * (-(t118 / 0.4e1 - t16 / 0.4e1) * t127 * t131 + (t106 + t86 * t114 + (-t16 * t115 + t48 * t152 + ((t115 * t25) - t126) * qJ(1)) * t26) * t20) - ((t103 * t3 + t138 * t4) * t47 + 0.2e1 * (t4 * (-(t49 * qJ(1)) / 0.2e1 + t48 * t10 * t79 - (0.2e1 * rSges(3,1) * t52 + rSges(3,3) * t103) * t123 / 0.2e1 + t14 * t128 / 0.2e1) + t3 * (-t52 * qJ(1) * t123 + t99 * t128 + (rSges(3,1) * t103 + (rSges(3,3) * t138)) * t123 / 0.2e1 + t16 * t102)) * t26) * qJD(1) / 0.2e1 + (-t14 * t4 + t16 * t3 + (-t25 * t3 + t4 * t60) * qJ(1)) * t151) * t58) * m(3) + (Icges(3,2) + Icges(2,3)) * ((t49 * t80) - t48 * t77 / 0.2e1 - t105 * t132 + (t85 * t114 + t20 * t133) * t26) + (-0.4e1 * t52 * t144 - 0.2e1 * (t142 * t32 + t36 * t143) * t27 * t46) * t53 * t125;];
tauc = t2(:);
